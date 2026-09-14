#include "loggerwidget.h"
#include <QTextCursor>
#include <QFileDialog>
#include <QCoreApplication>

// ============================================================================
// Logger Implementation
// ============================================================================

Logger& Logger::instance() {
    static Logger* instance = nullptr;
    if (!instance) {
        instance = new Logger();
    }
    return *instance;
}

Logger::Logger()
    : m_logFile(nullptr)
    , m_logStream(nullptr)
    , m_widget(nullptr)
    , m_initialized(false)
{
    // Don't do heavy initialization here
}

void Logger::ensureInitialized() {
    if (m_initialized) return;
    m_initialized = true;

    // Create logs directory
    QString logDir = QCoreApplication::applicationDirPath() + "/logs";
    QDir().mkpath(logDir);

    // Create log file
    QString timestamp = QDateTime::currentDateTime().toString("yyyy-MM-dd_HH-mm-ss");
    QString logPath = logDir + "/app_log_" + timestamp + ".txt";

    m_logFile = new QFile(logPath);
    if (m_logFile->open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Text)) {
        m_logStream = new QTextStream(m_logFile);
    }

    log(LogLevel::INFO, "=== Logger Initialized ===");
}

Logger::~Logger() {
    if (m_initialized) {
        log(LogLevel::INFO, "=== Logger Shutdown ===");
    }

    // Clean up timers
    QMutableMapIterator<QString, QElapsedTimer*> it(m_timers);
    while (it.hasNext()) {
        it.next();
        delete it.value();
    }
    m_timers.clear();

    // Clean up file
    if (m_logStream) {
        m_logStream->flush();
        delete m_logStream;
        m_logStream = nullptr;
    }
    if (m_logFile) {
        m_logFile->close();
        delete m_logFile;
        m_logFile = nullptr;
    }
}

void Logger::setWidget(LoggerWidget* widget) {
    QMutexLocker locker(&m_mutex);

    LoggerWidget* oldWidget = m_widget;
    m_widget = widget;

    // If we have a new widget and old was null, send buffered logs
    if (widget && !oldWidget) {
        // Send all buffered logs to the new widget
        for (const QString& log : m_logBuffer) {
            QMetaObject::invokeMethod(widget, "appendLog",
                                      Qt::QueuedConnection, Q_ARG(QString, log));
        }
    }
}

QStringList Logger::getBufferedLogs() {
    QMutexLocker locker(&m_mutex);
    return m_logBuffer;
}

QString Logger::getTimestamp() const {
    return QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss.zzz");
}

QString Logger::levelToString(LogLevel level) const {
    switch (level) {
    case LogLevel::DEBUG: return "DEBUG";
    case LogLevel::INFO:  return "INFO ";
    case LogLevel::WARN:  return "WARN ";
    case LogLevel::ERROR_LEVEL: return "ERROR";
    case LogLevel::TIMER: return "TIMER";
    default: return "UNKNOWN";
    }
}

void Logger::log(LogLevel level, const QString& message) {
    // Ensure logger is initialized (first call might trigger this)
    if (!m_initialized) {
        ensureInitialized();
    }

    QMutexLocker locker(&m_mutex);

    QString timestamp = getTimestamp();
    QString levelStr = levelToString(level);
    QString formattedMessage = QString("[%1] [%2] %3")
                                   .arg(timestamp)
                                   .arg(levelStr)
                                   .arg(message);

    // Add to buffer (with size limit)
    m_logBuffer.append(formattedMessage);
    if (m_logBuffer.size() > MAX_BUFFER_SIZE) {
        m_logBuffer.removeFirst();  // Remove oldest entry
    }

    writeToFile(formattedMessage);
    updateWidget(formattedMessage);

    // Also output to debug console
    if (level == LogLevel::ERROR_LEVEL || level == LogLevel::WARN) {
        qWarning().noquote() << formattedMessage;
    } else {
        qDebug().noquote() << formattedMessage;
    }
}

void Logger::startTimer(const QString& methodName) {
    if (!m_initialized) {
        ensureInitialized();
    }

    QMutexLocker locker(&m_mutex);

    if (m_timers.contains(methodName)) {
        delete m_timers[methodName];
    }

    QElapsedTimer* timer = new QElapsedTimer();
    timer->start();
    m_timers[methodName] = timer;
}

void Logger::stopTimer(const QString& methodName) {
    if (!m_initialized) {
        ensureInitialized();
    }

    QMutexLocker locker(&m_mutex);

    if (m_timers.contains(methodName)) {
        QElapsedTimer* timer = m_timers[methodName];
        qint64 elapsed = timer->elapsed();

        QString message = QString("%1 completed in %2 ms")
                              .arg(methodName)
                              .arg(elapsed);

        // Temporarily unlock to call log
        locker.unlock();
        log(LogLevel::TIMER, message);
        locker.relock();

        delete timer;
        m_timers.remove(methodName);
    }
}

void Logger::writeToFile(const QString& formattedMessage) {
    if (m_logStream && m_logFile && m_logFile->isOpen()) {
        *m_logStream << formattedMessage << "\n";
        m_logStream->flush();
    }
}

void Logger::updateWidget(const QString& formattedMessage) {
    if (m_widget) {
        QMetaObject::invokeMethod(m_widget, "appendLog",
                                  Qt::QueuedConnection, Q_ARG(QString, formattedMessage));
    }
}

// ============================================================================
// LoggerWidget Implementation
// ============================================================================

LoggerWidget::LoggerWidget(QWidget* parent) : QWidget(parent) {
    setWindowTitle("Application Log");
    resize(800, 600);

    QVBoxLayout* layout = new QVBoxLayout(this);

    m_logText = new QTextEdit(this);
    m_logText->setReadOnly(true);
    m_logText->setFont(QFont("Courier New", 9));

    m_logText->setStyleSheet(
        "QTextEdit {"
        "   background-color: #1e1e1e;"
        "   color: #d4d4d4;"
        "   font-family: 'Courier New', monospace;"
        "   font-size: 11px;"
        "}"
        );

    QHBoxLayout* buttonLayout = new QHBoxLayout();

    m_clearButton = new QPushButton("Clear", this);
    m_saveButton = new QPushButton("Save Log", this);

    buttonLayout->addStretch();
    buttonLayout->addWidget(m_clearButton);
    buttonLayout->addWidget(m_saveButton);

    layout->addWidget(m_logText);
    layout->addLayout(buttonLayout);

    connect(m_clearButton, &QPushButton::clicked, this, &LoggerWidget::onClearLogs);
    connect(m_saveButton, &QPushButton::clicked, this, &LoggerWidget::onSaveLogs);

    // Load all buffered logs first
    loadBufferedLogs();

    // Then register as widget (no need to send logs again since we already loaded them)
    Logger::instance().setWidget(this);
}

LoggerWidget::~LoggerWidget() {
    Logger::instance().setWidget(nullptr);
}

void LoggerWidget::loadBufferedLogs() {
    QStringList bufferedLogs = Logger::instance().getBufferedLogs();

    if (!bufferedLogs.isEmpty()) {
        m_logText->clear();
        for (const QString& log : bufferedLogs) {
            m_logText->append(log);
        }

        // Auto-scroll to bottom
        QTextCursor cursor = m_logText->textCursor();
        cursor.movePosition(QTextCursor::End);
        m_logText->setTextCursor(cursor);
    }
}

void LoggerWidget::appendLog(const QString& message) {
    // Append the message
    m_logText->append(message);

    // Auto-scroll to bottom
    QTextCursor cursor = m_logText->textCursor();
    cursor.movePosition(QTextCursor::End);
    m_logText->setTextCursor(cursor);
}

void LoggerWidget::onClearLogs() {
    m_logText->clear();
    Logger::instance().log(LogLevel::INFO, "Log window cleared");
}

void LoggerWidget::onSaveLogs() {
    QString fileName = QFileDialog::getSaveFileName(this,
                                                    "Save Log",
                                                    QCoreApplication::applicationDirPath() + "/logs/log_export.txt",
                                                    "Text Files (*.txt);;All Files (*)");

    if (!fileName.isEmpty()) {
        QFile file(fileName);
        if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            QTextStream stream(&file);
            stream << m_logText->toPlainText();
            file.close();

            Logger::instance().log(LogLevel::INFO,
                                   QString("Log exported to: %1").arg(fileName));
        } else {
            Logger::instance().log(LogLevel::ERROR_LEVEL,
                                   QString("Failed to export log to: %1").arg(fileName));
        }
    }
}


void Logger::closeLogFile()
{
    if (m_logStream) {
        m_logStream->flush();
        delete m_logStream;
        m_logStream = nullptr;
    }
    if (m_logFile) {
        m_logFile->close();
        delete m_logFile;
        m_logFile = nullptr;
    }
}

void Logger::setFileLogging(const QString& filePath)
{
    QMutexLocker locker(&m_mutex);

    // Закрываем текущий файл, если он открыт
    closeLogFile();

    // Создаём директорию, если нужно
    QFileInfo fileInfo(filePath);
    QDir().mkpath(fileInfo.absolutePath());

    // Открываем новый файл
    m_logFile = new QFile(filePath);
    if (m_logFile->open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Text)) {
        m_logStream = new QTextStream(m_logFile);
        m_fileLoggingEnabled = true;
        m_currentLogPath = filePath;

        locker.unlock();
        log(LogLevel::INFO, QString("File logging enabled: %1").arg(filePath));
        locker.relock();
    } else {
        m_fileLoggingEnabled = false;
        QString error = QString("Failed to open log file: %1").arg(filePath);
        locker.unlock();
        log(LogLevel::ERROR_LEVEL, error);
        locker.relock();
    }
}

void Logger::disableFileLogging()
{
    QMutexLocker locker(&m_mutex);

    if (m_fileLoggingEnabled) {
        locker.unlock();
        log(LogLevel::INFO, "File logging disabled");
        locker.relock();

        closeLogFile();
        m_fileLoggingEnabled = false;
        m_currentLogPath.clear();
    }
}