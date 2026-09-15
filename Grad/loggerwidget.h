#ifndef LOGGERWIDGET_H
#define LOGGERWIDGET_H

#include <QWidget>
#include <QTextEdit>
#include <QPushButton>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QFile>
#include <QTextStream>
#include <QDateTime>
#include <QMutex>
#include <QMutexLocker>
#include <QElapsedTimer>
#include <QApplication>
#include <QDebug>
#include <QDir>
#include <QStringList>

// Forward declaration
class LoggerWidget;

// Log levels
enum class LogLevel {
    DEBUG,
    INFO,
    WARN,
    ERROR_LEVEL,
    TIMER
};

// Singleton Logger class
class Logger {
public:
    static Logger& instance();

    void log(LogLevel level, const QString& message);
    void startTimer(const QString& methodName);
    void stopTimer(const QString& methodName);
    void setWidget(LoggerWidget* widget);

    // ============ НОВЫЕ МЕТОДЫ ============
    void setFileLogging(const QString& filePath);
    void disableFileLogging();
    bool isFileLoggingEnabled() const { return m_fileLoggingEnabled; }
    // =====================================

    // Get all buffered logs (non-const because we need to lock the mutex)
    QStringList getBufferedLogs();

private:
    Logger();
    ~Logger();
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    void ensureInitialized();
    QString levelToString(LogLevel level) const;
    QString getTimestamp() const;
    void writeToFile(const QString& formattedMessage);
    void updateWidget(const QString& formattedMessage);

    // ============ НОВЫЙ МЕТОД ============
    void closeLogFile();
    // =====================================

    QFile* m_logFile;
    QTextStream* m_logStream;
    mutable QMutex m_mutex;
    LoggerWidget* m_widget;
    QMap<QString, QElapsedTimer*> m_timers;
    bool m_initialized;

    // ============ НОВЫЕ ПЕРЕМЕННЫЕ ============
    bool m_fileLoggingEnabled;
    QString m_currentLogPath;
    // =========================================

    // Buffer to store all logs
    QStringList m_logBuffer;
    static const int MAX_BUFFER_SIZE = 10000;  // Maximum number of log entries to keep in memory
};

// LoggerWidget class
class LoggerWidget : public QWidget {
    Q_OBJECT

public:
    explicit LoggerWidget(QWidget* parent = nullptr);
    ~LoggerWidget();

public slots:
    void appendLog(const QString& message);

private slots:
    void onClearLogs();
    void onSaveLogs();

private:
    void loadBufferedLogs();

    QTextEdit* m_logText;
    QPushButton* m_clearButton;
    QPushButton* m_saveButton;
};

// RAII Timer class for automatic timing
class ScopedTimer {
public:
    ScopedTimer(const QString& methodName) : m_methodName(methodName) {
        Logger::instance().startTimer(m_methodName);
    }
    ~ScopedTimer() {
        Logger::instance().stopTimer(m_methodName);
    }
private:
    QString m_methodName;
};

#endif // LOGGERWIDGET_H