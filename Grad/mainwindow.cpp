 #include "mainwindow.h"
#include "ui_mainwindow.h"
#include <omp.h>

QPoint mPos;
QVector<QPoint> selectedEdge;
QVector<QPoint> trueEdge;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    m_dataModel(new AppDataModel(this))
{
    ui->setupUi(this);

    // Initialize logger widget - will be created lazily
    m_loggerWidget = nullptr;

    // Initialize widgets
    imCalculator = new ImageCalculator();
    imageWidget = new ImageShowcaseWidget();
    profileImageWidget = new ImageShowcaseWidget();
    getDataDialog = new Dialog();
    GrWid = new GraphWidget();

    // Set up connections
    connect(m_dataModel, &AppDataModel::currentImageChanged,
            this, &MainWindow::on_imageUpdated);
    connect(imCalculator, &ImageCalculator::throw_imageCalculator,
            this, &MainWindow::catch_ImageCalculator);
    connect(imageWidget, &ImageShowcaseWidget::imageClicked,
            this, &MainWindow::on_showcaseWidget_clicked);

    // Initialize parameters
    m_dataModel->setBlueMax(255);
    m_dataModel->setRedMax(255);
    m_dataModel->setRadius(10);
    m_dataModel->setSigma(5);

    // Set window properties
    this->setFixedSize(550, 50);

    // Update parameters from model
    iRad = m_dataModel->radius();
    dSigma = m_dataModel->sigma();
    dRedMax = m_dataModel->redMax();
    dBlueMax = m_dataModel->blueMax();

    // demonstrateBlurEffect();
    // demonstrateRadialBlurEffect();
    // demonstrateBlurEffectForRadiographicObject();

    // benchmarkConvolution();
    // benchmarkRefinement();
    generate_wth();
}

void MainWindow::generate_wth()
{

}

void MainWindow::demonstrateBlurEffect()
{
    Logger::instance().log(LogLevel::INFO, "=== Starting Blur Effect Demonstration ===");

    // ========================================================================
    // 1. Создание тестового изображения с резкой границей (окружность радиусом 80)
    // ========================================================================
    const int size = 200;
    const int radius = 50;
    const int center = size / 2;  // 100, 100

    QImage sharpImage(size, size, QImage::Format_ARGB32);
    sharpImage.fill(Qt::black);

    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            int dx = i - center;
            int dy = j - center;
            double dist = sqrt(dx*dx + dy*dy);

            if(dist <= radius) {
                int val = 255;
                sharpImage.setPixel(i, j, qRgb(val, val, val));
            } else {
                // Снаружи - чёрный (0)
                sharpImage.setPixel(i, j, qRgb(0, 0, 0));
            }
        }
    }

    // ========================================================================
    // 2. Размытие с разными параметрами сигмы
    // ========================================================================
    QVector<double> sigmaValues = {5.0, 10.0, 15.0};
    QVector<QImage> blurredImages;
    QVector<QString> labels;

    // Сохраняем исходное изображение
    blurredImages.push_back(sharpImage);
    labels.push_back("Оригинальное (σ=0)");

    for(double sigma : sigmaValues) {
        int kernelSize = 2 * int(4 * sigma + 0.5);
        if(kernelSize % 2 == 0) kernelSize++;  // Делаем нечётным

        Matrix2D<double> imageMat = ImageProcessor::fromGrayImage(sharpImage);
        Matrix2D<double> gaussKernel = ImageProcessor::getGauss(kernelSize, kernelSize, sigma);
        Matrix2D<double> blurredMat = ImageProcessor::convMat(imageMat, gaussKernel);
        blurredImages.push_back(ImageProcessor::toGrayImage(blurredMat));
        labels.push_back(QString("Размытое (σ=%1)").arg(sigma));
    }

    // ========================================================================
    // 3. Создание составного изображения с сеткой 1×4
    // ========================================================================
    const int thumbnailSize = 200;
    const int spacing = 10;
    const int totalWidth = thumbnailSize * 4 + spacing * 3;
    const int totalHeight = thumbnailSize + 50;  // +50 для подписей

    QImage compositeImage(totalWidth, totalHeight, QImage::Format_ARGB32);
    compositeImage.fill(Qt::white);

    QPainter painter(&compositeImage);
    painter.setPen(Qt::black);
    painter.setFont(QFont("Arial", 10));

    for(int i = 0; i < blurredImages.size(); i++) {
        int x = i * (thumbnailSize + spacing);

        // Масштабируем изображение до 200×200
        QImage scaled = blurredImages[i].scaled(thumbnailSize, thumbnailSize,
                                                Qt::KeepAspectRatio,
                                                Qt::SmoothTransformation);

        for(int i = 0; i < scaled.width(); i++) {
            for(int j = 0; j < scaled.height(); j++) {
                if(j == 200) scaled.setPixel({i,j},qRgb(0,255,255));
            }
        }
        // Рисуем изображение
        painter.drawImage(x, 0, scaled);

        // Рисуем рамку
        painter.drawRect(x, 0, thumbnailSize, thumbnailSize);

        // Рисуем подпись
        painter.drawText(x, thumbnailSize + 20, thumbnailSize, 25,
                         Qt::AlignCenter, labels[i]);
    }

    painter.end();

    // ========================================================================
    // 4. Отображение составного изображения
    // ========================================================================
    ImageShowcaseWidget* compositeWidget = new ImageShowcaseWidget();
    compositeWidget->setImage(compositeImage);
    compositeWidget->setWindowTitle("Сравнение размытия окружности с резкой границей (σ = 0, 5, 10, 15)");
    compositeWidget->resize(totalWidth + 20, totalHeight + 20);
    compositeWidget->show();

    // ========================================================================
    // 5. Построение горизонтальных профилей через центр
    // ========================================================================
    // Преобразуем все изображения в матрицы для построения профилей
    QVector<Matrix2D<double>> imageMatrices;
    imageMatrices.push_back(ImageProcessor::fromGrayImage(sharpImage));

    for(double sigma : sigmaValues) {
        int kernelSize = 2 * int(4 * sigma + 0.5);
        if(kernelSize % 2 == 0) kernelSize++;
        Matrix2D<double> imageMat = ImageProcessor::fromGrayImage(sharpImage);
        Matrix2D<double> gaussKernel = ImageProcessor::getGauss(kernelSize, kernelSize, sigma);
        Matrix2D<double> blurredMat = ImageProcessor::convMat(imageMat, gaussKernel);
        imageMatrices.push_back(blurredMat);
    }

    QVector<QString> profileLabels = {"Original (σ=0)", "σ=5", "σ=10", "σ=15"};
    QVector<QVector<QPointF>> horizontalProfiles;

    for(int idx = 0; idx < imageMatrices.size(); idx++) {
        const auto& image = imageMatrices[idx];
        QVector<QPointF> profile;

        // Горизонтальный профиль через центр (строка center = 100)
        for(int j = 0; j < size; j++) {
            double intensity = image[center][j];
            profile.push_back(QPointF(j, intensity));
        }
        horizontalProfiles.push_back(profile);
    }

    // Создаём окно с горизонтальными профилями
    GraphWidget* horizontalGraph = new GraphWidget();
    horizontalGraph->plotMultipleProfiles(horizontalProfiles, profileLabels);
    horizontalGraph->setAxisLabels("X coordinate (pixels)", "Intensity (0-255)");
    horizontalGraph->setTitle("Горизонтальные профили через центр изображения (y=100)");
    horizontalGraph->resize(800, 500);
    horizontalGraph->show();

    Logger::instance().log(LogLevel::INFO, "=== Blur Effect Demonstration Complete ===");
}
void MainWindow::demonstrateRadialBlurEffect()
{
    Logger::instance().log(LogLevel::INFO, "=== Starting Radial Blur Effect Demonstration (Cylinder with Cavity) ===");

    // ========================================================================
    // 1. Создание радиографического изображения (цилиндр с полостью)
    // ========================================================================
    const int size = 300;           // Размер изображения 300×300
    const int outerRadius = 100;    // Внешний радиус цилиндра
    const int innerRadius = 30;     // Радиус внутренней полости

    QImage radialImage(size, size, QImage::Format_ARGB32);
    radialImage.fill(Qt::black);

    // Максимальная интенсивность (соответствует минимальной толщине/плотности)
    const double maxIntensity = 255.0;

    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            int dx = i - 150;
            int dy = j - 150;
            int dx1 = i - 110;
            int dy1 = j - 110;
            double dist = sqrt(dx*dx + dy*dy);
            double dist1 = sqrt(dx1*dx1 + dy1*dy1);

            if(dist <= outerRadius) {
                double thickness = sqrt(outerRadius * outerRadius - dist * dist);

                // Если точка внутри полости, толщина уменьшается
                if(dist1 <= innerRadius) {
                    // Внутри полости - меньшая толщина (светлая область)
                    double cavityEffect = sqrt(innerRadius * innerRadius - dist1 * dist1);
                    thickness = thickness - cavityEffect*3;
                }

                // Нормируем на максимальную толщину (в центре сплошного цилиндра)
                double maxThickness = outerRadius;
                double intensity = (thickness / maxThickness) * maxIntensity;
                intensity = std::max(0.0, std::min(maxIntensity, intensity));

                int val = static_cast<int>(intensity);
                radialImage.setPixel(i, j, qRgb(val, val, val));
            } else {
                // Фон - чёрный (вне объекта)
                radialImage.setPixel(i, j, qRgb(0, 0, 0));
            }
        }
    }

    // ========================================================================
    // 2. Размытие с разными параметрами сигмы
    // ========================================================================
    QVector<double> sigmaValues = {3.0, 6.0, 9.0};
    QVector<QImage> blurredImages;
    QVector<QString> labels;

    // Сохраняем исходное изображение
    blurredImages.push_back(radialImage);
    labels.push_back("Оригинал (σ=0)");

    for(double sigma : sigmaValues) {
        int kernelSize = 2 * int(4 * sigma + 0.5);
        if(kernelSize % 2 == 0) kernelSize++;

        Matrix2D<double> imageMat = ImageProcessor::fromGrayImage(radialImage);
        Matrix2D<double> gaussKernel = ImageProcessor::getGauss(kernelSize, kernelSize, sigma);
        Matrix2D<double> blurredMat = ImageProcessor::convMat(imageMat, gaussKernel);
        blurredImages.push_back(ImageProcessor::toGrayImage(blurredMat));
        labels.push_back(QString("Размытие (σ=%1)").arg(sigma));
    }

    // ========================================================================
    // 3. Создание составного изображения с сеткой 1×4
    // ========================================================================
    const int thumbnailSize = 250;
    const int spacing = 15;
    const int totalWidth = thumbnailSize * 4 + spacing * 3;
    const int totalHeight = thumbnailSize + 50;

    QImage compositeImage(totalWidth, totalHeight, QImage::Format_ARGB32);
    compositeImage.fill(Qt::white);

    QPainter painter(&compositeImage);
    painter.setPen(Qt::black);
    painter.setFont(QFont("Arial", 10));

    for(int i = 0; i < blurredImages.size(); i++) {
        int x = i * (thumbnailSize + spacing);

        QImage scaled = blurredImages[i].scaled(thumbnailSize, thumbnailSize,
                                                Qt::KeepAspectRatio,
                                                Qt::SmoothTransformation);

        painter.drawImage(x, 0, scaled);
        painter.drawRect(x, 0, thumbnailSize, thumbnailSize);
        painter.drawText(x, thumbnailSize + 20, thumbnailSize, 25,
                         Qt::AlignCenter, labels[i]);
    }

    painter.end();

    // Отображаем составное изображение
    ImageShowcaseWidget* compositeWidget = new ImageShowcaseWidget();
    compositeWidget->setImage(compositeImage);
    compositeWidget->setWindowTitle("Радиографический объект: цилиндр с полостью (σ = 0, 3, 6, 9)");
    compositeWidget->resize(totalWidth + 20, totalHeight + 20);
    compositeWidget->show();

    // ========================================================================
    // 4. Построение профилей СУММИРОВАНИЕМ ПО СТОЛБЦАМ (лучевые суммы)
    //    Это имитирует реальную радиографию: каждый пиксель = интеграл по вертикали
    // ========================================================================

    // Преобразуем все изображения в матрицы
    QVector<Matrix2D<double>> imageMatrices;
    imageMatrices.push_back(ImageProcessor::fromGrayImage(radialImage));

    for(double sigma : sigmaValues) {
        int kernelSize = 2 * int(4 * sigma + 0.5);
        if(kernelSize % 2 == 0) kernelSize++;
        Matrix2D<double> imageMat = ImageProcessor::fromGrayImage(radialImage);
        Matrix2D<double> gaussKernel = ImageProcessor::getGauss(kernelSize, kernelSize, sigma);
        Matrix2D<double> blurredMat = ImageProcessor::convMat(imageMat, gaussKernel);
        imageMatrices.push_back(blurredMat);
    }

    // ========================================================================
    // 5. Горизонтальные профили (суммирование по Y для имитации лучевых сумм)
    //    Для каждого X суммируем все Y (интегрируем по вертикали)
    // ========================================================================
    QVector<QString> profileLabels = {"Original (σ=0)", "σ=3", "σ=6", "σ=9"};
    QVector<QVector<QPointF>> horizontalSumProfiles;

    for(int idx = 0; idx < imageMatrices.size(); idx++) {
        const auto& image = imageMatrices[idx];
        QVector<QPointF> profile;

        // Для каждого столбца X суммируем все значения по Y
        for(int x = 0; x < size; x++) {
            double columnSum = 0.0;
            for(int y = 0; y < size; y++) {
                columnSum += image[y][x];
            }
            // Нормируем для удобства отображения
            double normalizedSum = columnSum / (size * 0.5);
            profile.push_back(QPointF(x, normalizedSum));
        }
        horizontalSumProfiles.push_back(profile);
    }

    // Создаём окно с горизонтальными профилями (лучевые суммы)
    GraphWidget* horizontalSumGraph = new GraphWidget();
    horizontalSumGraph->plotMultipleProfiles(horizontalSumProfiles, profileLabels);
    horizontalSumGraph->setAxisLabels("X coordinate (pixels)", "Integrated intensity (ray sum)");
    horizontalSumGraph->setTitle("Горизонтальные профили (суммирование по Y) - имитация лучевых сумм");
    horizontalSumGraph->resize(900, 500);
    horizontalSumGraph->show();

    Logger::instance().log(LogLevel::INFO, "=== Radial Blur Effect Demonstration Complete ===");
}
void MainWindow::demonstrateBlurEffectForRadiographicObject()
{
    Logger::instance().log(LogLevel::INFO, "=== Starting Radiographic Object Blur Effect Demonstration ===");

    // ========================================================================
    // 1. Создание тестового изображения радиографического объекта (окружность с плавным изменением интенсивности)
    // ========================================================================
    const int size = 200;
    const int radius = 70;
    const int center = size / 2;  // 100, 100

    QImage radiographicImage(size, size, QImage::Format_ARGB32);
    radiographicImage.fill(Qt::black);

    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            int dx = i - center;
            int dy = j - center;
            double dist = sqrt(dx*dx + dy*dy);

            if(dist <= radius) {
                // Внутри окружности - интенсивность убывает от центра к краю
                double intensity = 255.0 * (1.0 - dist / radius);
                intensity = std::max(0.0, std::min(255.0, intensity));
                int val = static_cast<int>(intensity);
                radiographicImage.setPixel(i, j, qRgb(val, val, val));
            } else {
                radiographicImage.setPixel(i, j, qRgb(0, 0, 0));
            }
        }
    }

    // ========================================================================
    // 2. Размытие с разными параметрами сигмы
    // ========================================================================
    QVector<double> sigmaValues = {0.0, 3.0, 6.0, 9.0};
    QVector<Matrix2D<double>> blurredMatrices;
    QVector<QImage> blurredImages;

    for(double sigma : sigmaValues) {
        Matrix2D<double> imageMat = ImageProcessor::fromGrayImage(radiographicImage);

        if (sigma > 0) {
            int kernelSize = 2 * int(4 * sigma + 0.5);
            if(kernelSize % 2 == 0) kernelSize++;
            Matrix2D<double> gaussKernel = ImageProcessor::getGauss(kernelSize, kernelSize, sigma);
            Matrix2D<double> blurredMat = ImageProcessor::convMat(imageMat, gaussKernel);
            blurredMatrices.push_back(blurredMat);
            blurredImages.push_back(ImageProcessor::toGrayImage(blurredMat));
        } else {
            blurredMatrices.push_back(imageMat);
            blurredImages.push_back(radiographicImage);
        }
    }

    // ========================================================================
    // 3. Применение Лапласиана к каждому изображению
    // ========================================================================
    QVector<Matrix2D<double>> laplacianMatrices;
    QVector<QImage> laplacianImages;

    for(int idx = 0; idx < blurredMatrices.size(); idx++) {
        int kernelSize = 2 * int(4 * 2.0 + 0.5);
        if(kernelSize % 2 == 0) kernelSize++;
        Matrix2D<double> laplKernel = ImageProcessor::getLapl(kernelSize, kernelSize, 2.0);
        Matrix2D<double> laplacianMat = ImageProcessor::convMat(blurredMatrices[idx], laplKernel);
        laplacianMatrices.push_back(laplacianMat);
        laplacianImages.push_back(ImageProcessor::toBlueRedImage(laplacianMat, 255.0, 255.0));
    }

    // ========================================================================
    // 4. Создание составного изображения с сеткой 1×4 (с линией профиля)
    // ========================================================================
    const int thumbnailSize = 200;
    const int spacing = 10;
    const int totalWidth = thumbnailSize * 4 + spacing * 3;
    const int totalHeight = thumbnailSize + 80;

    QImage compositeImage(totalWidth, totalHeight, QImage::Format_ARGB32);
    compositeImage.fill(Qt::white);

    QPainter painter(&compositeImage);
    painter.setPen(Qt::black);
    painter.setFont(QFont("Arial", 9));

    QStringList labels = {"Оригинал (σ=0)", "Размытие (σ=3)", "Размытие (σ=6)", "Размытие (σ=9)"};

    // Координаты для линии профиля (горизонтальная линия через центр)
    int profileY = thumbnailSize / 2;  // 100 пикселей

    for(int i = 0; i < blurredImages.size(); i++) {
        int x = i * (thumbnailSize + spacing);

        QImage scaled = blurredImages[i].scaled(thumbnailSize, thumbnailSize,
                                                Qt::KeepAspectRatio,
                                                Qt::SmoothTransformation);

        // Рисуем изображение
        painter.drawImage(x, 0, scaled);

        // Рисуем красную горизонтальную линию там, где проходит профиль
        painter.setPen(QPen(Qt::yellow, 2, Qt::SolidLine));
        painter.drawLine(x, profileY, x + thumbnailSize - 1, profileY);

        // Рисуем рамку
        painter.setPen(QPen(Qt::black, 1));
        painter.drawRect(x, 0, thumbnailSize, thumbnailSize);

        // Рисуем подпись
        painter.setPen(Qt::black);
        painter.drawText(x, thumbnailSize + 5, thumbnailSize, 20,
                         Qt::AlignCenter, labels[i]);

        // Информация о сигме
        painter.drawText(x, thumbnailSize + 25, thumbnailSize, 20,
                         Qt::AlignCenter, QString("σ=%1").arg(sigmaValues[i]));
    }

    painter.end();

    ImageShowcaseWidget* compositeWidget = new ImageShowcaseWidget();
    compositeWidget->setImage(compositeImage);
    compositeWidget->setWindowTitle("Радиографический объект: исходное и размытые изображения (красная линия - профиль)");
    compositeWidget->resize(totalWidth + 20, totalHeight + 20);
    compositeWidget->show();

    // ========================================================================
    // 5. Лапласианы с линией профиля
    // ========================================================================
    QImage compositeLaplacian(totalWidth, totalHeight, QImage::Format_ARGB32);
    compositeLaplacian.fill(Qt::white);

    QPainter painterLapl(&compositeLaplacian);
    painterLapl.setPen(Qt::black);
    painterLapl.setFont(QFont("Arial", 9));

    QStringList laplLabels = {"Laplacian (σ=0)", "Laplacian (σ=3)", "Laplacian (σ=6)", "Laplacian (σ=9)"};

    for(int i = 0; i < laplacianImages.size(); i++) {
        int x = i * (thumbnailSize + spacing);

        QImage scaled = laplacianImages[i].scaled(thumbnailSize, thumbnailSize,
                                                  Qt::KeepAspectRatio,
                                                  Qt::SmoothTransformation);

        painterLapl.drawImage(x, 0, scaled);

        // Рисуем красную горизонтальную линию там, где проходит профиль
        painterLapl.setPen(QPen(Qt::yellow, 2, Qt::SolidLine));
        painterLapl.drawLine(x, profileY, x + thumbnailSize - 1, profileY);

        // Рисуем рамку
        painterLapl.setPen(QPen(Qt::black, 1));
        painterLapl.drawRect(x, 0, thumbnailSize, thumbnailSize);

        // Рисуем подпись
        painterLapl.setPen(Qt::black);
        painterLapl.drawText(x, thumbnailSize + 5, thumbnailSize, 20,
                             Qt::AlignCenter, laplLabels[i]);
        painterLapl.drawText(x, thumbnailSize + 25, thumbnailSize, 20,
                             Qt::AlignCenter, "Красный=положительные\n, Синий=отрицательные");
    }

    painterLapl.end();

    ImageShowcaseWidget* compositeLaplacianWidget = new ImageShowcaseWidget();
    compositeLaplacianWidget->setImage(compositeLaplacian);
    compositeLaplacianWidget->setWindowTitle("Результат применения Лапласиана");
    compositeLaplacianWidget->resize(totalWidth + 20, totalHeight + 20);
    compositeLaplacianWidget->show();

    // ========================================================================
    // 6. Построение исходных профилей (без Лапласиана)
    // ========================================================================
    QVector<QVector<QPointF>> originalProfiles;
    QVector<QString> profileLabels = {"Оригинал (σ=0)", "Размытие (σ=3)", "Размытие (σ=6)", "Размытие (σ=9)"};

    for(int idx = 0; idx < blurredMatrices.size(); idx++) {
        const auto& mat = blurredMatrices[idx];
        QVector<QPointF> profile;

        // Горизонтальный профиль через центр
        for(int j = 0; j < size; j++) {
            double intensity = mat[center][j];
            profile.push_back(QPointF(j, intensity));
        }
        originalProfiles.push_back(profile);
    }

    GraphWidget* originalProfileGraph = new GraphWidget();
    originalProfileGraph->plotMultipleProfiles(originalProfiles, profileLabels);
    originalProfileGraph->setAxisLabels("X coordinate (pixels)", "Intensity (0-255)");
    originalProfileGraph->setTitle("Исходные профили интенсивности через центр изображения (y=100)");
    originalProfileGraph->resize(800, 500);
    originalProfileGraph->show();

    // ========================================================================
    // 7. Построение профилей Лапласиана и поиск zero-crossing
    // ========================================================================
    QVector<QVector<QPointF>> laplacianProfiles;
    QVector<QVector<double>> zeroCrossings;

    for(int idx = 0; idx < laplacianMatrices.size(); idx++) {
        const auto& laplMat = laplacianMatrices[idx];
        QVector<QPointF> profile;
        QVector<double> zeros;

        // Горизонтальный профиль через центр
        for(int j = 0; j < size; j++) {
            double value = laplMat[center][j];
            profile.push_back(QPointF(j, value));
        }

        // Поиск zero-crossing
        for(int j = 0; j < size - 1; j++) {
            double v1 = laplMat[center][j];
            double v2 = laplMat[center][j + 1];

            if (v1 * v2 < 0) {
                double t = -v1 / (v2 - v1);
                double zeroX = j + t;
                zeros.push_back(zeroX);
            }
        }

        laplacianProfiles.push_back(profile);
        zeroCrossings.push_back(zeros);
    }

    GraphWidget* laplacianGraph = new GraphWidget();
    laplacianGraph->plotMultipleProfilesWithZeroCrossings(laplacianProfiles, profileLabels, zeroCrossings);
    laplacianGraph->setAxisLabels("X coordinate (pixels)", "Laplacian response");
    laplacianGraph->setTitle("Профили Лапласиана через центр изображения с отмеченными zero-crossing (красный пунктир)");
    laplacianGraph->resize(900, 600);
    laplacianGraph->show();

    Logger::instance().log(LogLevel::INFO, "=== Radiographic Object Blur Effect Demonstration Complete ===");
}


void MainWindow::benchmarkConvolution()
{
    Logger::instance().log(LogLevel::INFO, "=== Benchmark: Convolution ===");

    // Загружаем тестовое изображение
    QImage testImage = ImageProcessor::sampleTwoHollows();
    Matrix2D<double> imageMat = ImageProcessor::fromGrayImage(testImage);

    // Параметры для свертки (ядро Гаусса, сигма = 5, размер ядра 41x41)
    double sigma = 5.0;
    int kernelSize = 2 * int(4 * sigma + 0.5);
    if (kernelSize % 2 == 0) kernelSize++;
    Matrix2D<double> gaussKernel = ImageProcessor::getGauss(kernelSize, kernelSize, sigma);

    // Количество потоков для тестирования
    QVector<int> threadCounts = {1, 2, 4, 8};
    QVector<long long> times;
    const int iterations = 10;  // Количество повторений для усреднения

    Logger::instance().log(LogLevel::INFO, "| Threads | Time (ms) | Speedup | Efficiency |");
    Logger::instance().log(LogLevel::INFO, "|---------|-----------|---------|------------|");

    long long baseTime = 0;

    for (int threads : threadCounts) {
        omp_set_num_threads(threads);

        QElapsedTimer timer;
        timer.start();

        // Выполняем свертку iterations раз для усреднения
        for (int iter = 0; iter < iterations; iter++) {
            Matrix2D<double> result = ImageProcessor::convMat(imageMat, gaussKernel);
        }

        long long elapsed = timer.elapsed() / iterations;
        times.append(elapsed);

        if (threads == 1) {
            baseTime = elapsed;
        }

        double speedup = static_cast<double>(baseTime) / elapsed;
        double efficiency = (speedup / threads) * 100.0;

        Logger::instance().log(LogLevel::INFO,
                               QString("| %1       | %2         | %3     | %4%%        |")
                                   .arg(threads)
                                   .arg(elapsed)
                                   .arg(speedup, 0, 'f', 2)
                                   .arg(efficiency, 0, 'f', 1));
    }

    Logger::instance().log(LogLevel::INFO, "=== Benchmark Convolution Complete ===");
}

void MainWindow::benchmarkRefinement()
{
    Logger::instance().log(LogLevel::INFO, "=== Benchmark: PCR Refinement ===");

    // Загружаем тестовое изображение
    QImage testImage = ImageProcessor::sampleTwoHollows();
    m_dataModel->setCurrentImage(testImage);
    m_dataModel->setOriginalImage(testImage);

    int NX = testImage.width();
    int NY = testImage.height();

    // Параметры для ПЦР (как в test_002)
    double sigma1 = 4.0;
    double sigma2 = 8.0;
    double sigma0 = 8.0;
    double sigma01 = sqrt(sigma0 * sigma0 + 4.0 * 4.0);
    int NS = 2 * int(4 * sigma01 + 0.5);

    // Подготовка изображений
    Matrix2D<double> A0 = ImageProcessor::fromGrayImage(testImage);
    Matrix2D<double> B1 = ImageProcessor::convMat(A0, ImageProcessor::getGauss(NS, NS, sigma0));
    Matrix2D<double> B01 = ImageProcessor::convMat(A0, ImageProcessor::getGauss(NS, NS, sigma01));
    Matrix2D<double> A = B1;

    // Обновляем NS для max sigma
    double sigma_max = std::max(sigma1, sigma2);
    NS = 2 * int(4 * sigma_max + 0.5);

    // Вычисляем градиенты для получения направления
    Matrix2D<double> GxW1 = ImageProcessor::convMat(A, ImageProcessor::getXGradCore(NS, NS, sigma1));
    Matrix2D<double> GyW1 = ImageProcessor::convMat(A, ImageProcessor::getYGradCore(NS, NS, sigma1));

    int n0 = 78;
    int m0 = 83;

    double ex = GxW1[n0][m0];
    double ey = GyW1[n0][m0];
    double gradMag = sqrt(ex * ex + ey * ey);
    if (gradMag > 1e-6) {
        ex /= gradMag;
        ey /= gradMag;
    } else {
        ex = 1.0;
        ey = 0.0;
    }

    // Параметры для refineSinglePoint
    RefinementParameters params;
    params.A = A;
    params.B01 = B01;
    params.NX = NX;
    params.NY = NY;
    params.ex = ex;
    params.ey = ey;
    params.sigma0 = sigma0;
    params.sigma01 = sigma01;
    params.sigma1 = sigma1;
    params.sigma2 = sigma2;
    params.n_sigma = 10;
    params.n_myu = 200;
    params.NN = 100;
    params.otstup = 15;

    // Количество потоков для тестирования
    QVector<int> threadCounts = {1, 2, 4, 8};
    QVector<long long> times;

    Logger::instance().log(LogLevel::INFO, "| Threads | Time (ms) | Speedup | Efficiency |");
    Logger::instance().log(LogLevel::INFO, "|---------|-----------|---------|------------|");

    int iterations = 10;
    long long baseTime = 0;

    for (int threads : threadCounts) {
        omp_set_num_threads(threads);

        QElapsedTimer timer;
        timer.start();

        for (int i = 0; i < iterations; i++) {
            RefinementResult result = ImageProcessor::refineSinglePoint001(n0, m0, params);
        }

        long long elapsed = timer.elapsed() / iterations;
        times.append(elapsed);

        if (threads == 1) {
            baseTime = elapsed;
        }

        double speedup = static_cast<double>(baseTime) / elapsed;
        double efficiency = (speedup / threads) * 100.0;

        Logger::instance().log(LogLevel::INFO,
                               QString("| %1       | %2         | %3     | %4%%        |")
                                   .arg(threads)
                                   .arg(elapsed)
                                   .arg(speedup, 0, 'f', 2)
                                   .arg(efficiency, 0, 'f', 2));
    }

    Logger::instance().log(LogLevel::INFO, "=== Benchmark Refinement Complete ===");
}

void MainWindow::closeEvent(QCloseEvent *event)
{
    QString logFile = QCoreApplication::applicationDirPath() + "/logs/last_session.log";
    Logger::instance().setFileLogging(logFile);
    Logger::instance().log(LogLevel::INFO, "===Application shutdown===");

    // Небольшая задержка, чтобы лог записался
    QThread::msleep(100);

    event->accept();
}

MainWindow::~MainWindow()
{
    Logger::instance().log(LogLevel::INFO, "MainWindow destroyed");
    delete ui;
    delete profileImageWidget;
}

// ============================================================================
// File Operations
// ============================================================================

void MainWindow::on_actionopen_file_triggered()
{
    Logger::instance().log(LogLevel::INFO, "Open file dialog triggered");

    QFileDialog* openImDialog = new QFileDialog(this);
    openImDialog->setFileMode(QFileDialog::AnyFile);
    openImDialog->setNameFilter(tr("Images (*.png *.xpm *.jpg *.bmp *.tif *.tiff)"));

    QString fileName = openImDialog->getOpenFileName();
    if (fileName.isEmpty()) {
        Logger::instance().log(LogLevel::INFO, "File open cancelled by user");
        return;
    }

    Logger::instance().startTimer("Load Image");
    QImage newIm(fileName);
    Logger::instance().stopTimer("Load Image");

    if (newIm.isNull()) {
        Logger::instance().log(LogLevel::ERROR_LEVEL, QString("Failed to load image: %1").arg(fileName));
        QMessageBox::warning(this, "Error", "Failed to load image: " + fileName);
        return;
    }

    Logger::instance().log(LogLevel::INFO, QString("Image loaded: %1 (%2x%3)").arg(fileName).arg(newIm.width()).arg(newIm.height()));

    m_dataModel->setCurrentImage(newIm);
    m_dataModel->setOriginalImage(newIm);
    edgeSelectionMode = false;
    m_profileBuildingMode = false;
}

void MainWindow::on_actionsave_file_triggered()
{
    Logger::instance().log(LogLevel::INFO, "Save file dialog triggered");

    if (curImage.isNull()) {
        Logger::instance().log(LogLevel::ERROR_LEVEL, "Save failed: No image loaded");
        QMessageBox::warning(this, "Error", "No image to save");
        return;
    }

    QFileDialog* saveImDialog = new QFileDialog(this);
    saveImDialog->setFileMode(QFileDialog::AnyFile);
    saveImDialog->setNameFilter(tr("Images (*.png *.xpm *.jpg *.bmp *.tif *.tiff)"));
    QString saveFileName = saveImDialog->getSaveFileName();

    if (!saveFileName.isEmpty()) {
        Logger::instance().startTimer("Save Image");
        bool saved = curImage.save(saveFileName);
        Logger::instance().stopTimer("Save Image");

        if (saved) {
            Logger::instance().log(LogLevel::INFO, QString("Image saved: %1").arg(saveFileName));
        } else {
            Logger::instance().log(LogLevel::ERROR_LEVEL, QString("Failed to save image: %1").arg(saveFileName));
        }
    } else {
        Logger::instance().log(LogLevel::INFO, "Save cancelled by user");
    }
}

// ============================================================================
// Image Processing Operations
// ============================================================================

void MainWindow::on_actiongaussian_blur_triggered()
{
    Logger::instance().log(LogLevel::INFO, "Gaussian blur triggered");

    if (curImage.isNull()) {
        Logger::instance().log(LogLevel::ERROR_LEVEL, "Gaussian blur failed: No image loaded");
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    // Ask for sigma
    getDataDialog->setWindowTitle("Gaussian Blur - Sigma");
    getDataDialog->setPlaceholderText("Sigma value (e.g., 5.0)");
    getDataDialog->setDefaultValue(m_dataModel->sigma());

    if (!getDataDialog->exec()) {
        Logger::instance().log(LogLevel::INFO, "Gaussian blur cancelled by user");
        return;
    }
    double sigma = getDataDialog->getValue();

    if (sigma <= 0) {
        Logger::instance().log(LogLevel::ERROR_LEVEL, "Invalid sigma value");
        QMessageBox::warning(this, "Error", "Sigma must be positive");
        return;
    }

    // Ask for radius
    int defaultRadius = 8 * sigma;
    getDataDialog->setWindowTitle("Gaussian Blur - Radius");
    getDataDialog->setPlaceholderText(QString("Kernel radius (e.g., %1)").arg(defaultRadius));
    getDataDialog->setDefaultValue(defaultRadius);

    if (!getDataDialog->exec()) {
        Logger::instance().log(LogLevel::INFO, "Gaussian blur cancelled by user");
        return;
    }
    int radius = static_cast<int>(getDataDialog->getValue());

    if (radius <= 0 || radius % 2 == 0) {
        if (radius % 2 == 0) radius++; // Make odd
        Logger::instance().log(LogLevel::WARN,
                               QString("Radius adjusted to odd value: %1").arg(radius));
    }

    // Apply
    m_dataModel->setSigma(sigma);
    m_dataModel->setRadius(radius);
    iRad = m_dataModel->radius();
    dSigma = m_dataModel->sigma();

    Logger::instance().startTimer("Gaussian Blur");
    curImage = ImageProcessor::convImage(curImage,
                                         ImageProcessor::getGauss(iRad, iRad, dSigma));
    Logger::instance().stopTimer("Gaussian Blur");

    updateImageDisplay();
    im1 = curImage;

    Logger::instance().log(LogLevel::INFO,
                           QString("Gaussian blur applied (sigma=%1, radius=%2)").arg(dSigma).arg(iRad));
}

void MainWindow::on_actiongaussian_edge_detection_triggered()
{
    Logger::instance().log(LogLevel::INFO, "Gaussian edge detection triggered");

    if (curImage.isNull()) {
        Logger::instance().log(LogLevel::ERROR_LEVEL, "Edge detection failed: No image loaded");
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    getDataDialog->setWindowTitle("Write sigma coefficient value");
    getDataDialog->setPlaceholderText("Sigma value");
    getDataDialog->setDefaultValue(5.0);

    if (getDataDialog->exec()) {
        double sigma = getDataDialog->getValue();
        m_dataModel->setSigma(sigma);
        m_dataModel->setRadius(3 * m_dataModel->sigma());

        Logger::instance().startTimer("Gaussian Edge Detection");
        curImage = ImageProcessor::gaussianEdgeDetection(
            ImageProcessor::fromGrayImage(curImage),
            m_dataModel->sigma(),
            m_dataModel->radius(),
            attMat);
        Logger::instance().stopTimer("Gaussian Edge Detection");

        updateImageDisplay();
        edgeSelectionMode = true;

        Logger::instance().log(LogLevel::INFO, QString("Gaussian edge detection applied (sigma=%1)").arg(sigma));
    } else {
        Logger::instance().log(LogLevel::INFO, "Edge detection cancelled by user");
    }
}

void MainWindow::on_actionLaplacian_edge_detection_triggered()
{
    Logger::instance().log(LogLevel::INFO, "Laplacian edge detection triggered");

    if (curImage.isNull()) {
        Logger::instance().log(LogLevel::ERROR_LEVEL, "Laplacian failed: No image loaded");
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    getDataDialog->setWindowTitle("Write sigma coefficient value");
    getDataDialog->setPlaceholderText("Sigma value");
    getDataDialog->setDefaultValue(5.0);

    if (getDataDialog->exec()) {
        double sigma = getDataDialog->getValue();
        m_dataModel->setSigma(sigma);
        m_dataModel->setRadius(3 * m_dataModel->sigma());

        Logger::instance().startTimer("Laplacian Edge Detection");
        Matrix2D<double> mat1 = ImageProcessor::fromGrayImage(curImage);
        Matrix2D<double> mat2 = ImageProcessor::convMat(mat1,
                                                        ImageProcessor::getLapl(m_dataModel->radius(), m_dataModel->radius(), m_dataModel->sigma()));

        mat1 = ImageProcessor::findEdges(mat2, 3, attMat);
        curImage = ImageProcessor::toGrayImage(mat1);
        Logger::instance().stopTimer("Laplacian Edge Detection");

        updateImageDisplay();
        edgeSelectionMode = true;

        Logger::instance().log(LogLevel::INFO, QString("Laplacian edge detection applied (sigma=%1)").arg(sigma));
    } else {
        Logger::instance().log(LogLevel::INFO, "Laplacian cancelled by user");
    }
}

void MainWindow::on_actionsharpen_triggered()
{
    Logger::instance().log(LogLevel::INFO, "Sharpen filter triggered");

    if (curImage.isNull()) {
        Logger::instance().log(LogLevel::ERROR_LEVEL, "Sharpen failed: No image loaded");
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    Logger::instance().startTimer("Sharpen Filter");
    Matrix2D<double> kernel = {{-1, -1, -1},
                               {-1,  9, -1},
                               {-1, -1, -1}};
    curImage = ImageProcessor::convImage(curImage, kernel);
    Logger::instance().stopTimer("Sharpen Filter");

    updateImageDisplay();
    im1 = curImage;

    Logger::instance().log(LogLevel::INFO, "Sharpen filter applied");
}

void MainWindow::on_actiongradient_X_and_Y_triggered()
{
    Logger::instance().log(LogLevel::INFO, "Gradient X and Y triggered");

    if (curImage.isNull()) {
        Logger::instance().log(LogLevel::ERROR_LEVEL, "Gradient failed: No image loaded");
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    Logger::instance().startTimer("Gradient X and Y");

    Logger::instance().startTimer("Gradient X");
    Matrix2D<double> xGrad = ImageProcessor::fromGrayImage(curImage);
    xGrad = ImageProcessor::convMat(xGrad, ImageProcessor::getXGradCore(iRad, iRad, dSigma));
    ImageProcessor::toBlueRedImage(xGrad, 255., 255.).save("XGradImage.png");
    Logger::instance().stopTimer("Gradient X");

    Logger::instance().startTimer("Gradient Y");
    Matrix2D<double> yGrad = ImageProcessor::fromGrayImage(curImage);
    yGrad = ImageProcessor::convMat(yGrad, ImageProcessor::getYGradCore(iRad, iRad, dSigma));
    ImageProcessor::toBlueRedImage(yGrad, 255., 255.).save("YGradImage.png");
    Logger::instance().stopTimer("Gradient Y");

    Logger::instance().stopTimer("Gradient X and Y");

    Logger::instance().log(LogLevel::INFO, "Gradient images saved to XGradImage.png and YGradImage.png");
}

// ============================================================================
// Sample Images
// ============================================================================

void MainWindow::on_actionTwo_hollows_triggered()
{
    Logger::instance().log(LogLevel::INFO, "Sample 'Two Hollows' loaded");
    Logger::instance().startTimer("Generate Two Hollows");

    QImage setImage = ImageProcessor::sampleTwoHollows();
    m_dataModel->setCurrentImage(setImage);
    m_dataModel->setOriginalImage(setImage);
    curImage = setImage;
    edgeSelectionMode = false;
    m_profileBuildingMode = false;


    // First circle (bottom)
    int centerX1 = 100, centerY1 = 135, radius1 = 40;
    // Second circle (top)
    int centerX2 = 100, centerY2 = 65, radius2 = 30;

    trueEdge.clear();
    for(int i = 0; i < 200; i++){
        for(int j = 0; j < 200; j++){
            int dx1 = i - centerX1;
            int dy1 = j - centerY1;
            double dist1 = sqrt(static_cast<double>(dx1*dx1 + dy1*dy1));

            double dmin1 = dist1 - 0.7071;
            double dmax1 = dist1 + 0.7071;

            if(dmin1 <= radius1 && radius1 <= dmax1){
                trueEdge.push_back({i,j});
            }

            int dx2 = i - centerX2;
            int dy2 = j - centerY2;
            double dist2 = sqrt(static_cast<double>(dx2*dx2 + dy2*dy2));

            double dmin2 = dist2 - 0.7071;
            double dmax2 = dist2 + 0.7071;

            if(dmin2 <= radius2 && radius2 <= dmax2){
                trueEdge.push_back({i,j});
            }
        }
    }

    on_imageUpdated(setImage);
    Logger::instance().stopTimer("Generate Two Hollows");
    Logger::instance().log(LogLevel::INFO, QString("True edge points found: %1").arg(trueEdge.size()));
}

void MainWindow::on_actionTwo_hollows_big_triggered()
{
    Logger::instance().log(LogLevel::INFO, "Sample 'Two Hollows Big' loaded");
    Logger::instance().startTimer("Generate Two Hollows Big");

    QImage setImage = ImageProcessor::sampleTwoHollowsBig();
    m_dataModel->setCurrentImage(setImage);
    m_dataModel->setOriginalImage(setImage);
    edgeSelectionMode = false;
    m_profileBuildingMode = false;


    // First circle (bottom)
    int centerX1 = 200, centerY1 = 200, radius1 = 100;
    // Second circle (top)
    int centerX2 = 200, centerY2 = 375, radius2 = 75;

    trueEdge.clear();
    for(int i = 0; i < 400; i++){
        for(int j = 0; j < 550; j++){
            int dx1 = i - centerX1;
            int dy1 = j - centerY1;
            double dist1 = sqrt(static_cast<double>(dx1*dx1 + dy1*dy1));

            double dmin1 = dist1 - 0.7071;
            double dmax1 = dist1 + 0.7071;

            if(dmin1 <= radius1 && radius1 <= dmax1){
                trueEdge.push_back({i,j});
            }

            int dx2 = i - centerX2;
            int dy2 = j - centerY2;
            double dist2 = sqrt(static_cast<double>(dx2*dx2 + dy2*dy2));

            double dmin2 = dist2 - 0.7071;
            double dmax2 = dist2 + 0.7071;

            if(dmin2 <= radius2 && radius2 <= dmax2){
                trueEdge.push_back({i,j});
            }
        }
    }

    on_imageUpdated(setImage);
    Logger::instance().stopTimer("Generate Two Hollows Big");
    Logger::instance().log(LogLevel::INFO, QString("True edge points found: %1").arg(trueEdge.size()));
}

// ============================================================================
// Profile Operations
// ============================================================================

void MainWindow::on_actiondraw_profile_triggered()
{
    Logger::instance().log(LogLevel::INFO, "Draw profile triggered");

    if (curImage.isNull()) {
        Logger::instance().log(LogLevel::ERROR_LEVEL, "Profile failed: No image loaded");
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    Logger::instance().startTimer("Draw Vertical Profile");
    QVector<double> x, y;
    int centerX = curImage.width() / 2;

    for (int i = 0; i < curImage.height(); i++) {
        x.push_back(i);
        y.push_back(qGray(curImage.pixel(centerX, i)));
    }

    GrWid->plotGraph(x, y);
    GrWid->setAxisLabels("Y coordinate", "Intensity");
    GrWid->setTitle(QString("Вертикальный Профиль в x = %1").arg(centerX));
    GrWid->show();
    Logger::instance().stopTimer("Draw Vertical Profile");

    Logger::instance().log(LogLevel::INFO, QString("Vertical profile drawn at x=%1").arg(centerX));
}

void MainWindow::on_actionProfileBetweenPoints_triggered()
{
    Logger::instance().log(LogLevel::INFO, "Profile between points triggered");

    if (curImage.isNull()) {
        Logger::instance().log(LogLevel::ERROR_LEVEL, "Profile failed: No image loaded");
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    if (profilePoints.size() < 2) {
        m_profileBuildingMode = true;
        profilePoints.clear();
        Logger::instance().log(LogLevel::INFO, "Profile mode activated - click two points");
    } else {
        Logger::instance().startTimer("Build Profile Between Points");

        QVector<double> profile = ImageProcessor::buildProfileBetweenPoints(
            profilePoints[0], profilePoints[1], curImage);

        if (profile.isEmpty()) {
            Logger::instance().log(LogLevel::ERROR_LEVEL, "Failed to build profile between points");
            QMessageBox::warning(this, "Error", "Failed to build profile");
            Logger::instance().stopTimer("Build Profile Between Points");
            return;
        }

        // Рисуем линию профиля на изображении
        QImage imageWithLine = drawProfileLineOnImage(curImage, profilePoints[0], profilePoints[1], profile);

        // Показываем изображение с линией в отдельном окне
        ImageShowcaseWidget* profileImageWidget = new ImageShowcaseWidget();
        profileImageWidget->setImage(imageWithLine);
        profileImageWidget->setWindowTitle("Линия профиля на изображении");
        profileImageWidget->resize(400, 400);
        profileImageWidget->show();

        // Строим график профиля
        QVector<double> x(profile.size());
        for (int i = 0; i < profile.size(); i++) {
            x[i] = i;
        }

        GrWid->plotGraph(x, profile);
        GrWid->setAxisLabels("Distance along line (pixels)", "Intensity");
        GrWid->setTitle(QString("Profile from (%1,%2) to (%3,%4)")
                            .arg(profilePoints[0].x()).arg(profilePoints[0].y())
                            .arg(profilePoints[1].x()).arg(profilePoints[1].y()));
        GrWid->show();

        double minVal = *std::min_element(profile.begin(), profile.end());
        double maxVal = *std::max_element(profile.begin(), profile.end());
        double sum = 0.0;
        for (double v : profile) sum += v;
        double mean = sum / profile.size();

        Logger::instance().stopTimer("Build Profile Between Points");
        Logger::instance().log(LogLevel::INFO, QString("Profile stats - Length: %1, Min: %2, Max: %3, Mean: %4")
                                                   .arg(profile.size()).arg(minVal).arg(maxVal).arg(mean));

        profilePoints.clear();
        m_profileBuildingMode = false;
    }
}

void MainWindow::on_actionClearProfilePoints_triggered()
{
    Logger::instance().log(LogLevel::INFO, "Profile points cleared");
    profilePoints.clear();
    m_profileBuildingMode = false;
}

void MainWindow::onImageClickedForProfile(const QPoint& imagePosition)
{
    if (!m_profileBuildingMode) return;

    profilePoints.append(imagePosition);

    if (profilePoints.size() == 1) {
        Logger::instance().log(LogLevel::INFO, QString("First profile point: (%1,%2)").arg(imagePosition.x()).arg(imagePosition.y()));
    } else if (profilePoints.size() == 2) {
        Logger::instance().startTimer("Build Profile");

        QVector<double> profile = ImageProcessor::buildProfileBetweenPoints(
            profilePoints[0], profilePoints[1], curImage);

        if (!profile.isEmpty()) {
            // Рисуем линию профиля на изображении
            QImage imageWithLine = drawProfileLineOnImage(curImage, profilePoints[0], profilePoints[1], profile);

            // Показываем изображение с линией в отдельном окне
            ImageShowcaseWidget* profileImageWidget = new ImageShowcaseWidget();
            profileImageWidget->setImage(imageWithLine);
            profileImageWidget->setWindowTitle("Линия профиля на изображении");
            profileImageWidget->resize(400, 400);
            profileImageWidget->show();

            // Строим график профиля
            QVector<double> x(profile.size());
            for (int i = 0; i < profile.size(); i++) x[i] = i;

            GrWid->plotGraph(x, profile);
            GrWid->setAxisLabels("Distance along line (pixels)", "Intensity");
            GrWid->setTitle(QString("Profile from (%1,%2) to (%3,%4)")
                                .arg(profilePoints[0].x()).arg(profilePoints[0].y())
                                .arg(profilePoints[1].x()).arg(profilePoints[1].y()));
            GrWid->show();
        }

        Logger::instance().stopTimer("Build Profile");
        Logger::instance().log(LogLevel::INFO, "Profile between points displayed");

        profilePoints.clear();
        m_profileBuildingMode = false;
    }
}

// ============================================================================
// Statistics Operations
// ============================================================================

void MainWindow::on_actionShowStatistics_triggered()
{
    Logger::instance().log(LogLevel::INFO, "Show statistics triggered");

    if (curImage.isNull()) {
        Logger::instance().log(LogLevel::ERROR_LEVEL, "Statistics failed: No image loaded");
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    Logger::instance().startTimer("Compute Statistics");
    auto stats = ImageProcessor::computeImageStatistics(curImage, attMat);
    Logger::instance().stopTimer("Compute Statistics");

    Logger::instance().log(LogLevel::INFO, QString("Statistics: Mean=%1, Range=[%2, %3]").arg(stats.mean).arg(stats.min).arg(stats.max));
    showStatisticsDialog(stats);
}

void MainWindow::on_actionExportStatistics_triggered()
{
    Logger::instance().log(LogLevel::INFO, "Export statistics triggered");

    if (curImage.isNull()) {
        Logger::instance().log(LogLevel::ERROR_LEVEL, "Export failed: No image loaded");
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    QString fileName = QFileDialog::getSaveFileName(this, "Save Statistics",
                                                    "", "Text Files (*.txt);;All Files (*)");

    if (fileName.isEmpty()) {
        Logger::instance().log(LogLevel::INFO, "Export statistics cancelled");
        return;
    }

    Logger::instance().startTimer("Export Statistics");
    auto stats = ImageProcessor::computeImageStatistics(curImage, attMat);
    bool success = ImageProcessor::saveStatisticsToFile(stats, fileName);
    Logger::instance().stopTimer("Export Statistics");

    if (success) {
        Logger::instance().log(LogLevel::INFO, QString("Statistics exported to: %1").arg(fileName));
        QMessageBox::information(this, "Success", "Statistics saved to " + fileName);
    } else {
        Logger::instance().log(LogLevel::ERROR_LEVEL, QString("Failed to export statistics to: %1").arg(fileName));
        QMessageBox::warning(this, "Error", "Failed to save statistics");
    }
}

void MainWindow::on_actionShowLogWindow_triggered()
{
    if (!m_loggerWidget) {
        m_loggerWidget = new LoggerWidget();
    }
    m_loggerWidget->show();
    m_loggerWidget->raise();

    Logger::instance().log(LogLevel::INFO, "Log window opened");
}

// ============================================================================
// Edge Selection and Refinement
// ============================================================================

void MainWindow::on_showcaseWidget_clicked(const QPoint &imagePosition)
{
    mPos.setX(imagePosition.x());
    mPos.setY(imagePosition.y());

    // Handle profile building mode
    if (m_profileBuildingMode) {
        onImageClickedForProfile(imagePosition);
        return;
    }

    // Handle edge selection mode
    if (edgeSelectionMode && !attMat.empty()) {
        Logger::instance().startTimer("Edge Selection");
        resetEdgeSelection();

        selectedEdge.clear();
        for(int i = 0; i < attMat.size(); i++) {
            for(int j = 0; j < attMat[0].size(); j++) {
                if(attMat[i][j] == static_cast<int>(attribute::isEdge)) {
                    attMat[i][j] = static_cast<int>(attribute::isSelectedEdge);
                    selectedEdge.push_back({i,j}); }
            }
        }

        // ImageProcessor::selectEdge(mPos,attMat,selectedEdge);
        Logger::instance().log(LogLevel::INFO, QString("Edge selected at (%1,%2) - %3 edge points").arg(imagePosition.x()).arg(imagePosition.y()).arg(selectedEdge.size()));

        m_dataModel->setSelectedEdge(selectedEdge);

        // Highlight selected edge in cyan
        for (int i = 0; i < selectedEdge.size(); i++) {
            curImage.setPixel(selectedEdge[i], qRgb(0, 255, 255));
        }

        updateImageDisplay();
        Logger::instance().stopTimer("Edge Selection");
    }
}

void MainWindow::resetEdgeSelection()
{
    if (attMat.empty()) return;

    for (int i = 0; i < curImage.width(); i++) {
        for (int j = 0; j < curImage.height(); j++) {
            if (attMat[i][j] == static_cast<int>(attribute::isSelectedEdge)) {
                attMat[i][j] = static_cast<int>(attribute::isEdge);
                curImage.setPixel(i, j, qRgb(255, 255, 255));
            }
        }
    }
}

void MainWindow::on_actiontest_triggered()
{
    Logger::instance().log(LogLevel::INFO, "=== Edge Refinement Test (test_001) Started ===");
    Logger::instance().startTimer("Edge Refinement Test 001");

    if (!edgeSelectionMode) {
        Logger::instance().log(LogLevel::ERROR_LEVEL, "Test failed: Edge selection mode is off");
        QMessageBox::warning(this, "Error", "Edge selection mode is off. Run edge detection first.");
        Logger::instance().stopTimer("Edge Refinement Test 001");
        return;
    }

    // Collect all selected edge points
    QVector<QPoint> selectedEdgePoints = selectedEdge;

    if (selectedEdgePoints.isEmpty()) {
        Logger::instance().log(LogLevel::ERROR_LEVEL, "Test failed: No edge points selected");
        QMessageBox::warning(this, "Error", "No edge points selected");
        Logger::instance().stopTimer("Edge Refinement Test 001");
        return;
    }

    Logger::instance().log(LogLevel::INFO, QString("Refining %1 edge points").arg(selectedEdgePoints.size()));

    // Create images for visualization
    QImage trueEdgeImage = m_dataModel->originalImage();
    QImage gradientEdgeImage = m_dataModel->originalImage();
    QImage pcrEdgeImage = m_dataModel->originalImage();

    // Image dimensions
    int NX = m_dataModel->originalImage().width();
    int NY = m_dataModel->originalImage().height();

    // Gaussian blur parameters
    double sigma0 = 8.0;
    double sigma01 = sqrt(sigma0 * sigma0 + 4.0 * 4.0);
    int NS = 2 * int(4 * sigma01 + 0.5);

    // Two different sigma values for edge detection
    double sigma1 = 4.0;
    double sigma2 = 8.0;
    double sigma_max = max(sigma1, sigma2);

    // Prepare matrices
    Logger::instance().log(LogLevel::DEBUG, "Preparing matrices for refinement");
    Matrix2D<double> A0 = ImageProcessor::fromGrayImage(m_dataModel->originalImage());
    Matrix2D<double> B1 = ImageProcessor::convMat(A0, ImageProcessor::getGauss(NS, NS, sigma0));
    Matrix2D<double> B01 = ImageProcessor::convMat(A0, ImageProcessor::getGauss(NS, NS, sigma01));
    Matrix2D<double> A = B1;

    // Update kernel size for max sigma
    NS = 2 * int(4 * sigma_max + 0.5);

    // Calculate gradients
    Matrix2D<double> B = A;
    Matrix2D<double> GxW1 = ImageProcessor::convMat(B, ImageProcessor::getXGradCore(NS, NS, sigma1));
    Matrix2D<double> GyW1 = ImageProcessor::convMat(B, ImageProcessor::getYGradCore(NS, NS, sigma1));

    // Profile parameters
    int n_sigma = 10;
    int n_myu = 200;
    double sigma_myu = min(sigma1, sigma2);
    int NN = 100;
    int otstup = 15;

    RefinementParameters params;
    params.A = A;
    params.B01 = B01;
    params.NN = NN;
    params.NX = NX;
    params.NY = NY;
    params.n_myu = n_myu;
    params.n_sigma = n_sigma;
    params.otstup = otstup;
    params.sigma0 = sigma0;
    params.sigma01 = sigma01;
    params.sigma1 = sigma1;
    params.sigma2 = sigma2;

    // Thread-safe storage for results
    QVector<QPointF> refinedPositions;
    refinedPositions.resize(selectedEdgePoints.size());

    QVector<double> shiftsX, shiftsY;
    shiftsX.resize(selectedEdgePoints.size());
    shiftsY.resize(selectedEdgePoints.size());

#pragma omp for
    for (int pointIdx = 0; pointIdx < selectedEdgePoints.size(); pointIdx++) {
        int n0 = selectedEdgePoints[pointIdx].x();
        int m0 = selectedEdgePoints[pointIdx].y();

        double ex = GxW1[n0][m0];
        double ey = GyW1[n0][m0];
        double gradMag = sqrt(ex * ex + ey * ey);

        if (gradMag > 1e-6) {
            ex /= gradMag;
            ey /= gradMag;
        } else {
            ex = 1.0;
            ey = 0.0;
        }

        // Create local copy of params for each thread (thread-safe)
        RefinementParameters localParams = params;
        localParams.ex = ex;
        localParams.ey = ey;

        RefinementResult result = ImageProcessor::refineSinglePoint001(n0, m0, localParams);
        double n_new = result.refinedPosition.x();
        double m_new = result.refinedPosition.y();

        // Thread-safe assignment to specific index
        refinedPositions[pointIdx] = QPointF(n_new, m_new);
        shiftsX[pointIdx] = n_new - n0;
        shiftsY[pointIdx] = m_new - m0;

        if ((pointIdx + 1) % 10 == 0) {
            Logger::instance().log(LogLevel::DEBUG, QString("Refined %1/%2 points").arg(pointIdx + 1).arg(selectedEdgePoints.size()));
        }
    }

    // Draw true edge (green)
    for(int i = 0; i < trueEdge.size(); i++) {
        int x = trueEdge[i].x();
        int y = trueEdge[i].y();
        if (x >= 0 && x < trueEdgeImage.width() && y >= 0 && y < trueEdgeImage.height()) {
            trueEdgeImage.setPixel(x, y, qRgb(0, 255, 0));
        }
    }

    // Draw gradient edges (blue)
    for(int i = 0; i < selectedEdgePoints.size(); i++) {
        int x = selectedEdgePoints[i].x();
        int y = selectedEdgePoints[i].y();
        if (x >= 0 && x < gradientEdgeImage.width() && y >= 0 && y < gradientEdgeImage.height()) {
            gradientEdgeImage.setPixel(x, y, qRgb(0, 0, 255));
        }
    }

    // Draw PCR refined edges (red)
    for(int i = 0; i < refinedPositions.size(); i++) {
        int x = static_cast<int>(refinedPositions[i].x() + 0.5);
        int y = static_cast<int>(refinedPositions[i].y() + 0.5);
        if (x >= 0 && x < pcrEdgeImage.width() && y >= 0 && y < pcrEdgeImage.height()) {
            pcrEdgeImage.setPixel(x, y, qRgb(255, 0, 0));
        }
    }

    // Calculate statistics
    double sumShiftX = 0.0, sumShiftY = 0.0;
    double maxShiftX = 0.0, maxShiftY = 0.0;
    double minShiftX = 1e10, minShiftY = 1e10;

    for (int i = 0; i < selectedEdgePoints.size(); i++) {
        double absShiftX = std::abs(shiftsX[i]);
        double absShiftY = std::abs(shiftsY[i]);

        sumShiftX += absShiftX;
        sumShiftY += absShiftY;

        if (absShiftX > maxShiftX) maxShiftX = absShiftX;
        if (absShiftY > maxShiftY) maxShiftY = absShiftY;

        if (absShiftX < minShiftX) minShiftX = absShiftX;
        if (absShiftY < minShiftY) minShiftY = absShiftY;
    }

    double avgShiftX = (selectedEdgePoints.size() > 0) ? sumShiftX / selectedEdgePoints.size() : 0;
    double avgShiftY = (selectedEdgePoints.size() > 0) ? sumShiftY / selectedEdgePoints.size() : 0;

    // Show statistics
    Logger::instance().log(LogLevel::INFO, "=== REFINEMENT SUMMARY ===");
    Logger::instance().log(LogLevel::INFO, QString("Total points processed: %1").arg(refinedPositions.size()));
    Logger::instance().log(LogLevel::INFO, QString("Average absolute shift: ΔX = %1 pixels, ΔY = %2 pixels").arg(avgShiftX, 0, 'f', 3).arg(avgShiftY, 0, 'f', 3));
    Logger::instance().log(LogLevel::INFO, QString("Maximum absolute shift: ΔX = %1 pixels, ΔY = %2 pixels").arg(maxShiftX, 0, 'f', 3).arg(maxShiftY, 0, 'f', 3));
    Logger::instance().log(LogLevel::INFO, QString("Minimum absolute shift: ΔX = %1 pixels, ΔY = %2 pixels").arg(minShiftX, 0, 'f', 3).arg(minShiftY, 0, 'f', 3));

    // Calculate average shift magnitude
    double avgMagnitude = 0.0;
    double maxMagnitude = 0.0;
    double minMagnitude = 1e10;
    for (int i = 0; i < selectedEdgePoints.size(); i++) {
        double magnitude = std::sqrt(shiftsX[i] * shiftsX[i] + shiftsY[i] * shiftsY[i]);
        avgMagnitude += magnitude;
        if (magnitude > maxMagnitude) maxMagnitude = magnitude;
        if (magnitude < minMagnitude) minMagnitude = magnitude;
    }
    avgMagnitude /= selectedEdgePoints.size();

    Logger::instance().log(LogLevel::INFO, QString("Average shift magnitude: %1 pixels").arg(avgMagnitude, 0, 'f', 3));
    Logger::instance().log(LogLevel::INFO, QString("Maximum shift magnitude: %1 pixels").arg(maxMagnitude, 0, 'f', 3));
    Logger::instance().log(LogLevel::INFO, QString("Minimum shift magnitude: %1 pixels").arg(minMagnitude, 0, 'f', 3));

    // Show final images
    ImageShowcaseWidget* trueEdgeWidget = new ImageShowcaseWidget();
    trueEdgeWidget->setImage(trueEdgeImage);
    trueEdgeWidget->setWindowTitle("Истинная граница (зелёный)");
    trueEdgeWidget->resize(400, 400);
    trueEdgeWidget->show();

    ImageShowcaseWidget* gradientEdgeWidget = new ImageShowcaseWidget();
    gradientEdgeWidget->setImage(gradientEdgeImage);
    gradientEdgeWidget->setWindowTitle("Градиентная граница (синий)");
    gradientEdgeWidget->resize(400, 400);
    gradientEdgeWidget->show();

    ImageShowcaseWidget* pcrEdgeWidget = new ImageShowcaseWidget();
    pcrEdgeWidget->setImage(pcrEdgeImage);
    pcrEdgeWidget->setWindowTitle("Граница ПЦР (красный)");
    pcrEdgeWidget->resize(400, 400);
    pcrEdgeWidget->show();

    // Combined image with all three
    QImage combinedImage = m_dataModel->originalImage();

    for(int i = 0; i < trueEdge.size(); i++) {
        int x = trueEdge[i].x();
        int y = trueEdge[i].y();
        if (x >= 0 && x < combinedImage.width() && y >= 0 && y < combinedImage.height()) {
            combinedImage.setPixel(x, y, qRgb(0, 255, 0));
        }
    }

    for(int i = 0; i < selectedEdgePoints.size(); i++) {
        int x = selectedEdgePoints[i].x();
        int y = selectedEdgePoints[i].y();
        if (x >= 0 && x < combinedImage.width() && y >= 0 && y < combinedImage.height()) {
            // Don't overwrite green
            if (qGreen(combinedImage.pixel(x, y)) != 255) {
                combinedImage.setPixel(x, y, qRgb(0, 0, 255));
            }
        }
    }

    for(int i = 0; i < refinedPositions.size(); i++) {
        int x = static_cast<int>(refinedPositions[i].x() + 0.5);
        int y = static_cast<int>(refinedPositions[i].y() + 0.5);
        if (x >= 0 && x < combinedImage.width() && y >= 0 && y < combinedImage.height()) {
            combinedImage.setPixel(x, y, qRgb(255, 0, 0));
        }
    }

    ImageShowcaseWidget* combinedWidget = new ImageShowcaseWidget();
    combinedWidget->setImage(combinedImage);
    combinedWidget->setWindowTitle("Сравнение: зелёный — истина, синий — градиент, красный — ПЦР");
    combinedWidget->resize(400, 400);
    combinedWidget->show();

    Logger::instance().stopTimer("Edge Refinement Test 001");
    Logger::instance().log(LogLevel::INFO, "=== Edge Refinement Test Completed ===");
}

void MainWindow::on_actiontest_002_triggered()
{
    Logger::instance().log(LogLevel::INFO, "=== Edge Refinement Test (test_002) Started ===");
    Logger::instance().startTimer("Edge Refinement Test 002");

    // Generate test image
    QImage setImage = ImageProcessor::sampleTwoHollows();
    m_dataModel->setCurrentImage(setImage);
    m_dataModel->setOriginalImage(setImage);

    // Image dimensions
    int NX = setImage.width();
    int NY = setImage.height();

    // Gaussian blur parameters
    double sigma0 = 8.0;
    double sigma01 = sqrt(sigma0 * sigma0 + 4.0 * 4.0);
    int NS = 2 * int(4 * sigma01 + 0.5);

    // Two different sigma values for edge detection
    double sigma1 = 4.0;
    double sigma2 = 8.0;
    double sigma_max = max(sigma1, sigma2);

    Logger::instance().log(LogLevel::DEBUG, QString("Parameters: sigma0=%1, sigma1=%2, sigma2=%3, NS=%4").arg(sigma0).arg(sigma1).arg(sigma2).arg(NS));

    // Convert images to matrices and apply Gaussian blur
    Matrix2D<double> A0 = ImageProcessor::fromGrayImage(m_dataModel->currentImage());
    Matrix2D<double> B1 = ImageProcessor::convMat(A0, ImageProcessor::getGauss(NS, NS, sigma0));
    Matrix2D<double> B01 = ImageProcessor::convMat(A0, ImageProcessor::getGauss(NS, NS, sigma01));
    Matrix2D<double> A = B1;

    // Update kernel size for max sigma
    NS = 2 * int(4 * sigma_max + 0.5);

    // Perform Gaussian edge detection
    curImage = ImageProcessor::gaussianEdgeDetection(
        ImageProcessor::fromGrayImage(curImage), sigma_max, NS, attMat);
    edgeSelectionMode = true;
    updateImageDisplay();

    // Select test point for refinement
    int n0 = 78, m0 = 83;

    // Calculate gradients
    Matrix2D<double> B = A;
    Matrix2D<double> GxW1 = ImageProcessor::convMat(B, ImageProcessor::getXGradCore(NS, NS, sigma1));
    Matrix2D<double> GyW1 = ImageProcessor::convMat(B, ImageProcessor::getYGradCore(NS, NS, sigma1));

    double ex = GxW1[n0][m0];
    double ey = GyW1[n0][m0];
    double gradMag = sqrt(ex * ex + ey * ey);

    if (gradMag > 1e-6) {
        ex /= gradMag;
        ey /= gradMag;
    } else {
        ex = 1.0;
        ey = 0.0;
    }

    Logger::instance().log(LogLevel::DEBUG, QString("Gradient at (%1,%2): ex=%3, ey=%4, |grad|=%5").arg(n0).arg(m0).arg(ex).arg(ey).arg(gradMag));

    //*****************************************
    // PROFILE CONSTRUCTION
    //*****************************************

    QVector<double> prof1, prof2;
    int n_sigma = 10;
    int n_myu = 200;
    int x0 = n0, y0 = m0;
    double sigma_myu = min(sigma1, sigma2);
    QVector<double> mu;

    // Generate mu values
    for(int s = 0; s < n_myu; s++) {
        double val = -n_sigma*sigma_myu + s*(2*n_sigma*sigma_myu)/n_myu;
        mu.push_back(val);
    }

    // Build profiles along gradient direction
    for(int s = 0; s < n_myu; s++) {
        double prof1Val = 0.0;
        double prof2Val = 0.0;
        for(int n = 0; n < NX; n++) {
            for(int m = 0; m < NY; m++) {
                double mult = 0.0;
                double dx = (x0 - n + mu[s]*ex);
                double dy = (y0 - m + mu[s]*ey);
                double rad = dx*dx + dy*dy;
                mult = (rad/sigma1/sigma1 - 2.0) * exp(-1./2.*rad/sigma1/sigma1);
                prof1Val += A[n][m]*mult;
                prof2Val += B01[n][m]*mult;
            }
        }
        prof1.push_back(prof1Val);
        prof2.push_back(prof2Val);
    }

    // Find zero crossings
    App_Stats y2_stats;
    QVector<double> y2_for_stats;
    for(int i = 0; i < n_myu; i++) y2_for_stats.push_back(prof2[i]);
    y2_stats.gather_stats(y2_for_stats);
    int nmumax = max(y2_stats.x_max, y2_stats.x_min);
    int nmumin = min(y2_stats.x_max, y2_stats.x_min);
    int N0 = n_myu;

    int nL_Zero = N0/2, nR_Zero = N0/2;
    for(int n = nmumin; n < nmumax && n < n_myu - 1; n++) {
        if(prof1[n]*prof1[n+1] < 0) nL_Zero = n;
        if(prof2[n]*prof2[n+1] < 0) nR_Zero = n;
    }

    int N_Zero = (nL_Zero + nR_Zero)/2;

    // Find profile bounds by counting increases/decreases
    int decrCounter = 0, incrCounter = 0;
    int XL1 = 0, XR1 = 0;

    // Search right bound
    for(int i = N_Zero + 1; i < n_myu - 1; i++) {
        if(prof2[i] < prof2[i+1]) {
            incrCounter++;
            incrCounter = min(incrCounter, 5);
        }
        if(prof2[i] > prof2[i+1]) {
            decrCounter++;
            decrCounter = min(decrCounter, 5);
        }
        if(decrCounter == 5 && incrCounter == 5) {
            XR1 = min(i, n_myu - 35);
            break;
        }
    }
    if (XR1 == 0) XR1 = n_myu - 35;

    // Search left bound
    incrCounter = 0, decrCounter = 0;
    for(int i = N_Zero - 1; i > 1; i--) {
        if(prof2[i] < prof2[i+1]) {
            incrCounter++;
            incrCounter = min(incrCounter, 5);
        }
        if(prof2[i] > prof2[i+1]) {
            decrCounter++;
            decrCounter = min(decrCounter, 5);
        }
        if(decrCounter == 5 && incrCounter == 5) {
            XL1 = max(i, 20);
            break;
        }
    }
    if (XL1 == 0) XL1 = 20;

    Logger::instance().log(LogLevel::DEBUG, QString("Profile bounds: XL1=%1, XR1=%2").arg(XL1).arg(XR1));

    // Prepare first profile for interpolation
    QVector<QPointF> yP1, yP2;
    for(int i = XL1; i < XR1; i++) {
        yP1.push_back({ 1. * i, prof1[i]} );
    }

    int NN = 100;

    // Reinterpolate first profile
    QVector<QPointF> yyP1 = ImageProcessor::reInterpolateProfile(yP1, NN);
    QVector<double> yy1;
    for(int s1 = 0; s1 < NN; s1++) {
        yy1.push_back(yyP1[s1].y());
    }

    // Search for optimal second profile bounds
    int otstup = 15;
    double minraz = INT_MAX;
    int best_XL2 = XL1, best_XR2 = XR1;

    // Prepare full second profile for interpolation
    QVector<QPointF> yP2_full;
    for(int i = 0; i < n_myu; i++) {
        yP2_full.push_back({ 1. * i, prof2[i]} );
    }

    for(int XL2 = XL1; XL2 > XL1 - otstup; XL2--) {
        for(int XR2 = XR1; XR2 < XR1 + otstup; XR2++) {
            if (XL2 < 0 || XR2 >= n_myu || XL2 >= XR2) continue;

            QVector<QPointF> yyP2 = ImageProcessor::reInterpolateProfile(yP2_full, XL2, XR2, NN);
            QVector<double> yy2;
            for(int s1 = 0; s1 < NN; s1++) {
                yy2.push_back(yyP2[s1].y());
            }

            double m1 = 0.0, m2 = 0.0, D1 = 0.0, D2 = 0.0;
            for(int i = 0; i < NN; i++) {
                m1 += yy1[i];
                m2 += yy2[i];
            }
            m1 /= NN;
            m2 /= NN;

            for(int i = 0; i < NN; i++) {
                D1 += (yy1[i] - m1) * (yy1[i] - m1);
                D2 += (yy2[i] - m2) * (yy2[i] - m2);
            }
            D1 = sqrt(D1 / NN);
            D2 = sqrt(D2 / NN);

            if (D1 < 1e-10) D1 = 1.0;
            if (D2 < 1e-10) D2 = 1.0;

            double MINRAZ = 0.0;
            for(int i = 0; i < NN; i++) {
                double norm1 = (yy1[i] - m1) / D1;
                double norm2 = (yy2[i] - m2) / D2;
                double diff = norm1 - norm2;
                MINRAZ += diff * diff;
            }

            if(minraz > MINRAZ) {
                minraz = MINRAZ;
                best_XL2 = XL2;
                best_XR2 = XR2;
            }
        }
    }

    Logger::instance().log(LogLevel::DEBUG, QString("Best bounds: XL2=%1, XR2=%2, minraz=%3").arg(best_XL2).arg(best_XR2).arg(minraz));

    // Final processing with optimal bounds
    QVector<QPointF> yyP2_final = ImageProcessor::reInterpolateProfile(yP2_full, best_XL2, best_XR2, NN);
    QVector<double> yy2_final;
    for(int s1 = 0; s1 < NN; s1++) {
        yy2_final.push_back(yyP2_final[s1].y());
    }

    // Final normalization
    double m1_f = 0.0, m2_f = 0.0, D1_f = 0.0, D2_f = 0.0;
    for(int i = 0; i < NN; i++) {
        m1_f += yy1[i];
        m2_f += yy2_final[i];
    }
    m1_f /= NN;
    m2_f /= NN;

    for(int i = 0; i < NN; i++) {
        D1_f += (yy1[i] - m1_f) * (yy1[i] - m1_f);
        D2_f += (yy2_final[i] - m2_f) * (yy2_final[i] - m2_f);
    }
    D1_f = sqrt(D1_f / NN);
    D2_f = sqrt(D2_f / NN);

    if (D1_f < 1e-10) D1_f = 1.0;
    if (D2_f < 1e-10) D2_f = 1.0;

    // Normalized original profiles
    QVector<double> yk1, yk2;
    for(int s = 0; s < n_myu; s++) {
        yk1.push_back((prof1[s] - m1_f) / D1_f);
        yk2.push_back((prof2[s] - m2_f) / D2_f);
    }

    // Calculate shift coefficients
    double xa, ya, xb, yb;
    xa = XL1; xb = XR1;
    ya = 1. * best_XL2 - XL1;
    yb = 1. * best_XR2 - XR1;

    // Find intersection point
    double x00 = (xa*yb - xb*ya) / (yb - ya);
    double mu00 = -n_sigma*sigma_myu + x00*(2*n_sigma*sigma_myu)/n_myu;

    // Transform mu for second profile
    double FFF1_1 = mu00;      // Shift coefficient
    double FFF1_0 = (yb - ya) / (xb - xa);  // Scale coefficient

    QVector<double> mu_new;
    for(int i = 0; i < n_myu; i++) {
        double val = (mu[i] - FFF1_1) * (1 - FFF1_0) + FFF1_1;
        mu_new.push_back(val);
    }

    // Calculate new refined position
    double n_new = n0 + ex * FFF1_1;
    double m_new = m0 + ey * FFF1_1;

    Logger::instance().log(LogLevel::INFO, QString("Refinement result: (%1,%2) -> (%3,%4)").arg(n0).arg(m0).arg(n_new).arg(m_new));
    Logger::instance().log(LogLevel::INFO, QString("Shift: %1, Scale: %2").arg(FFF1_1).arg(FFF1_0));

    // ========================================================================
    // ПОСТРОЕНИЕ ВСЕХ ГРАФИКОВ В КОНЦЕ
    // ========================================================================

    // 1. График: Исходные профили Лапласиана (до сведения)
    GraphWidget* grwid2 = new GraphWidget();
    QVector<double> xIndices;
    for(int i = 0; i < n_myu; i++) xIndices.push_back(i);

    QCustomPlot* plot2 = grwid2->getPlot();
    plot2->addGraph();
    plot2->graph(0)->setData(xIndices, prof1);
    plot2->graph(0)->setPen(QPen(Qt::blue, 2));
    plot2->graph(0)->setName("Исходный профиль (σ=4)");

    plot2->addGraph();
    plot2->graph(1)->setData(xIndices, prof2);
    plot2->graph(1)->setPen(QPen(Qt::red, 2));
    plot2->graph(1)->setName("Размытый профиль (σ=8)");

    plot2->xAxis->setLabel("Номер точки (индекс)");
    plot2->yAxis->setLabel("Значение Лапласиана");
    plot2->legend->setVisible(true);
    plot2->rescaleAxes(true);
    plot2->replot();
    grwid2->setWindowTitle("Профили Лапласиана до сведения");
    grwid2->resize(800, 500);
    grwid2->show();

    // 2. График: Профили с отмеченными границами (XL1, XR1, best_XL2, best_XR2)
    GraphWidget* grwid_bounds = new GraphWidget();
    QCustomPlot* plot_bounds = grwid_bounds->getPlot();
    plot_bounds->addGraph();
    plot_bounds->graph(0)->setData(xIndices, prof1);
    plot_bounds->graph(0)->setPen(QPen(Qt::blue, 2));
    plot_bounds->graph(0)->setName("Исходный профиль (σ=4)");
    plot_bounds->graph(0)->setScatterStyle(QCPScatterStyle::ssNone);

    plot_bounds->addGraph();
    plot_bounds->graph(1)->setData(xIndices, prof2);
    plot_bounds->graph(1)->setPen(QPen(Qt::red, 2));
    plot_bounds->graph(1)->setName("Размытый профиль (σ=8)");
    plot_bounds->graph(1)->setScatterStyle(QCPScatterStyle::ssNone);

    // Точки границ первого профиля (XL1, XR1)
    QCPGraph* bounds1 = plot_bounds->addGraph();
    bounds1->setData(QVector<double>{double(XL1), double(XR1)},
                     QVector<double>{prof1[XL1], prof1[XR1]});
    bounds1->setPen(QPen(Qt::darkGreen, 2));
    bounds1->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 10));
    bounds1->setLineStyle(QCPGraph::lsNone);
    bounds1->setName("Пределы приграничного интервала  исходного профиля (XL1, XR1)");

    // Точки границ второго профиля (best_XL2, best_XR2)
    QCPGraph* bounds2 = plot_bounds->addGraph();
    bounds2->setData(QVector<double>{double(best_XL2), double(best_XR2)},
                     QVector<double>{prof2[best_XL2], prof2[best_XR2]});
    bounds2->setPen(QPen(Qt::darkRed, 2));
    bounds2->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDiamond, 10));
    bounds2->setLineStyle(QCPGraph::lsNone);
    bounds2->setName("Оптимальные пределы приграничного интервала размытого профиля (XL2, XR2)");

    // Точка zero-crossing первого профиля
    int firstBoundaryPoint = 100;
    for(int i = 1; i < prof1.size(); i++) {
        if(prof1[i-1]*prof1[i] <= 0) firstBoundaryPoint = i;
    }
    plot_bounds->addGraph();
    plot_bounds->graph(4)->setData(QVector<double>{double(firstBoundaryPoint)},
                                   QVector<double>{prof1[firstBoundaryPoint]});
    plot_bounds->graph(4)->setPen(QPen(Qt::cyan, 2));
    plot_bounds->graph(4)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCross, 12));
    plot_bounds->graph(4)->setName("Исходная позиция границы (zero-crossing)");

    // Точка zero-crossing второго профиля (найденная ПЦР) - в пространстве mu
    // Находим индекс в массиве mu для FFF1_1
    int secondBoundaryIndex = -1;
    for(int i = 0; i < n_myu - 1; i++) {
        if (mu[i] <= FFF1_1 && mu[i+1] >= FFF1_1) {
            secondBoundaryIndex = i;
            break;
        }
    }

    double secondBoundaryValue = 0;
    double secondBoundaryX = 0;
    if (secondBoundaryIndex != -1) {
        // Находим значение профиля в точке FFF1_1 для второго профиля (yk2)
        double t = (FFF1_1 - mu[secondBoundaryIndex]) / (mu[secondBoundaryIndex + 1] - mu[secondBoundaryIndex]);
        secondBoundaryValue = yk2[secondBoundaryIndex] * (1 - t) + yk2[secondBoundaryIndex + 1] * t;
        // Находим соответствующий индекс в исходном профиле (для отображения на графике)
        secondBoundaryX = (secondBoundaryIndex + t) / (double)n_myu * (xIndices.last() - xIndices.first()) + xIndices.first();
    }

    plot_bounds->addGraph();
    plot_bounds->graph(5)->setData(QVector<double>{secondBoundaryIndex}, QVector<double>{prof2[secondBoundaryIndex]});
    plot_bounds->graph(5)->setPen(QPen(Qt::magenta, 2));
    plot_bounds->graph(5)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssStar, 14));
    plot_bounds->graph(5)->setName(QString("Позиция границы ПЦР (μ = %1)").arg(FFF1_1, 0, 'f', 3));

    plot_bounds->xAxis->setLabel("Номер точки (индекс)");
    plot_bounds->yAxis->setLabel("Значение Лапласиана");
    plot_bounds->legend->setVisible(true);
    plot_bounds->rescaleAxes(true);
    plot_bounds->replot();
    grwid_bounds->setWindowTitle("Профили с найденными оптимальными границами");
    grwid_bounds->resize(800, 600);
    grwid_bounds->show();

    // 3. График: Окончательно сведённые профили (после применения сдвига и масштабирования)
    GraphWidget* grwid_aligned = new GraphWidget();
    QCustomPlot* plot_aligned = grwid_aligned->getPlot();

    plot_aligned->addGraph();
    plot_aligned->graph(0)->setData(mu, yk1);
    plot_aligned->graph(0)->setPen(QPen(Qt::blue, 2));
    plot_aligned->graph(0)->setName("Исходный профиль (σ=4)");

    plot_aligned->addGraph();
    plot_aligned->graph(1)->setData(mu_new, yk2);
    plot_aligned->graph(1)->setPen(QPen(Qt::red, 2));
    plot_aligned->graph(1)->setName("Размытый профиль после сжатия (σ=8)");

    // Отмечаем точку границы ПЦР на втором профиле
    double muValueAtShift = 0;
    for(int i = 0; i < n_myu - 1; i++) {
        if (mu[i] <= FFF1_1 && mu[i+1] >= FFF1_1) {
            double t = (FFF1_1 - mu[i]) / (mu[i+1] - mu[i]);
            muValueAtShift = yk2[i] * (1 - t) + yk2[i+1] * t;
            break;
        }
    }

    plot_aligned->addGraph();
    plot_aligned->graph(2)->setData(QVector<double>{FFF1_1}, QVector<double>{muValueAtShift});
    plot_aligned->graph(2)->setPen(QPen(Qt::magenta, 2));
    plot_aligned->graph(2)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssStar, 14));
    plot_aligned->graph(2)->setName(QString("Позиция границы ПЦР (μ = %1)").arg(FFF1_1, 0, 'f', 3));

    plot_aligned->xAxis->setLabel("Параметр μ (расстояние вдоль градиента)");
    plot_aligned->yAxis->setLabel("Нормализованное значение Лапласиана");
    plot_aligned->legend->setVisible(true);
    plot_aligned->rescaleAxes(true);
    plot_aligned->replot();
    grwid_aligned->setWindowTitle("Окончательно сведённые профили после применения сжатия");
    grwid_aligned->resize(800, 500);
    grwid_aligned->show();

    // Create visualization of original and refined points
    QImage tmp = m_dataModel->originalImage();
    QImage whiteIm = tmp;
    for(int i = 0; i < tmp.width(); i++){
        for(int j = 0; j < tmp.height(); j++){
            if(qGray(whiteIm.pixel(i,j)) != 0.) whiteIm.setPixel(i,j,qRgb(255,255,255));
            else whiteIm.setPixel(i,j,qRgb(0,0,0));
        }
    }

    // Mark original and refined points
    whiteIm.setPixel(n0, m0, qRgb(0, 0, 255));      // Blue - original point
    whiteIm.setPixel(static_cast<int>(n_new + 0.5), static_cast<int>(m_new + 0.5), qRgb(255, 0, 0));  // Red - refined point

    // Show final result
    ImageShowcaseWidget* resultWidget = new ImageShowcaseWidget();
    resultWidget->setImage(whiteIm);
    resultWidget->setWindowTitle("Результат уточнения положения границы (синий — исходная, красный — уточнённая)");
    resultWidget->resize(400, 400);
    resultWidget->show();

    Logger::instance().stopTimer("Edge Refinement Test 002");
    Logger::instance().log(LogLevel::INFO, "=== Edge Refinement Test 002 Completed ===");
}

// ============================================================================
// Image Calculator
// ============================================================================

void MainWindow::on_actionimage_calculator_triggered()
{
    Logger::instance().log(LogLevel::INFO, "Image calculator triggered");

    if (curImage.isNull()) {
        Logger::instance().log(LogLevel::ERROR_LEVEL, "Calculator failed: No image loaded");
        QMessageBox::warning(this, "Error", "No image loaded");
        return;
    }

    imCalculator->setImage1(curImage, "Current Image");
    imCalculator->show();
    imCalculator->setWindowTitle("Image Calculator");
}

void MainWindow::catch_ImageCalculator(const QImage& image1, const QImage& image2,
                                       QString operation, bool newWindow, bool floatResult)
{
    Q_UNUSED(floatResult);

    Logger::instance().log(LogLevel::INFO, QString("Image Calculator: %1 operation").arg(operation));
    Logger::instance().startTimer(QString("ImageCalc_%1").arg(operation));

    if (image1.isNull()) {
        Logger::instance().log(LogLevel::ERROR_LEVEL, "Image Calculator: Image 1 is null");
        return;
    }

    QImage res(image1.width(), image1.height(), QImage::Format_ARGB32);
    Matrix2D<double> ddMat1 = ImageProcessor::fromGrayImage(image1);
    Matrix2D<double> ddMat2 = ImageProcessor::fromGrayImage(image2);

    if (operation == "Add") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::AddWithMaxClamp<double>{}));
    } else if (operation == "Subtract") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::SubtractWithZeroClamp<double>{}));
    } else if (operation == "Multiply") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::Multiply<double>{}));
    } else if (operation == "Divide") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::Divide<double>{}));
    } else if (operation == "AND") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::BitwiseAND<double>{}));
    } else if (operation == "OR") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::BitwiseOR<double>{}));
    } else if (operation == "XOR") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::BitwiseXOR<double>{}));
    } else if (operation == "Min") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::Min<double>{}));
    } else if (operation == "Max") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::Max<double>{}));
    } else if (operation == "Average") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::Average<double>{}));
    } else if (operation == "Difference") {
        res = ImageProcessor::toGrayImage(
            ImageProcessor::elementWiseOperation(ddMat1, ddMat2, MatrixLambdas::Difference<double>{}));
    } else if (operation == "Copy") {
        res = image1;
    }

    Logger::instance().stopTimer(QString("ImageCalc_%1").arg(operation));

    if (newWindow) {
        ImageShowcaseWidget* showIm = new ImageShowcaseWidget();
        showIm->setAttribute(Qt::WA_DeleteOnClose);
        showIm->setImage(res);
        showIm->setWindowTitle("Image Calculator Result");
        showIm->show();
        Logger::instance().log(LogLevel::INFO, "Calculator result displayed in new window");
    } else {
        updateImageDisplay();
        Logger::instance().log(LogLevel::INFO, "Calculator result displayed in main window");
    }
}

// ============================================================================
// Image Display and Updates
// ============================================================================

void MainWindow::on_imageUpdated(const QImage &newImage)
{
    // Clear matrices
    clearMatrix2D(attMat);
    clearMatrix2D(matForTest);
    clearMatrix2D(imMat);

    curImage = newImage;
    im1 = curImage;
    imMat = ImageProcessor::fromGrayImage(curImage);

    int width = curImage.width();
    int height = curImage.height();

    // Resize matrices
    attMat.resize(width);
    for (int i = 0; i < width; i++) {
        attMat[i].resize(height, 0);
    }

    matForTest.resize(width);
    for (int i = 0; i < width; i++) {
        matForTest[i].resize(height, 0.0);
    }

    updateImageDisplay();

    Logger::instance().log(LogLevel::DEBUG, QString("Image updated: %1x%2").arg(width).arg(height));
}

void MainWindow::updateImageDisplay()
{
    imageWidget->setImage(curImage);
    imageWidget->show();
}

void MainWindow::showStatisticsDialog(const ImageProcessor::ImageStatistics& stats)
{
    QMessageBox msgBox;
    msgBox.setWindowTitle("Image Statistics");
    msgBox.setText(stats.toString());
    msgBox.setTextFormat(Qt::PlainText);
    msgBox.setStandardButtons(QMessageBox::Ok);

    // Add export buttons
    QPushButton* exportStatsBtn = msgBox.addButton("Export Statistics", QMessageBox::ActionRole);
    QPushButton* exportHistBtn = msgBox.addButton("Export Histogram", QMessageBox::ActionRole);

    msgBox.exec();

    if (msgBox.clickedButton() == exportStatsBtn) {
        on_actionExportStatistics_triggered();
    }
}
void MainWindow::showProfileDialog(const QVector<double>& profile, const QString& title)
{
    if (profile.isEmpty()) return;

    QDialog dialog(this);
    dialog.setWindowTitle(title);
    dialog.resize(600, 400);

    QVBoxLayout* layout = new QVBoxLayout(&dialog);

    GraphWidget* graph = new GraphWidget(&dialog);
    QVector<double> x(profile.size());
    for (int i = 0; i < profile.size(); i++) x[i] = i;
    graph->plotGraph(x, profile);
    graph->setAxisLabels("Distance (pixels)", "Intensity");
    graph->setTitle(title);

    layout->addWidget(graph);

    QPushButton* closeBtn = new QPushButton("Close", &dialog);
    layout->addWidget(closeBtn);

    connect(closeBtn, &QPushButton::clicked, &dialog, &QDialog::accept);

    dialog.exec();
}

// ============================================================================
// Unused Actions (Placeholders)
// ============================================================================

void MainWindow::on_actionpcr_edge_detection_triggered()
{
    Logger::instance().log(LogLevel::INFO, "=== Edge Refinement Test (test_001) Started ===");
    Logger::instance().startTimer("Edge Refinement Test 001");

    if (!edgeSelectionMode) {
        Logger::instance().log(LogLevel::ERROR_LEVEL, "Test failed: Edge selection mode is off");
        QMessageBox::warning(this, "Error", "Edge selection mode is off. Run edge detection first.");
        Logger::instance().stopTimer("Edge Refinement Test 001");
        return;
    }

    // Collect all selected edge points
    QVector<QPoint> selectedEdgePoints = selectedEdge;

    if (selectedEdgePoints.isEmpty()) {
        Logger::instance().log(LogLevel::ERROR_LEVEL, "Test failed: No edge points selected");
        QMessageBox::warning(this, "Error", "No edge points selected");
        Logger::instance().stopTimer("Edge Refinement Test 001");
        return;
    }

    Logger::instance().log(LogLevel::INFO, QString("Refining %1 edge points").arg(selectedEdgePoints.size()));

    // Create images for visualization
    QImage trueEdgeImage = m_dataModel->originalImage();
    QImage gradientEdgeImage = m_dataModel->originalImage();
    QImage pcrEdgeImage = m_dataModel->originalImage();

    // Image dimensions
    int NX = m_dataModel->originalImage().width();
    int NY = m_dataModel->originalImage().height();

    // Gaussian blur parameters
    double sigma0 = 8.0;
    double sigma01 = sqrt(sigma0 * sigma0 + 4.0 * 4.0);
    int NS = 2 * int(4 * sigma01 + 0.5);

    // Two different sigma values for edge detection
    double sigma1 = 4.0;
    double sigma2 = 8.0;
    double sigma_max = max(sigma1, sigma2);

    // Prepare matrices
    Logger::instance().log(LogLevel::DEBUG, "Preparing matrices for refinement");
    Matrix2D<double> A0 = ImageProcessor::fromGrayImage(m_dataModel->originalImage());
    Matrix2D<double> B1 = ImageProcessor::convMat(A0, ImageProcessor::getGauss(NS, NS, sigma0));
    Matrix2D<double> B01 = ImageProcessor::convMat(A0, ImageProcessor::getGauss(NS, NS, sigma01));
    Matrix2D<double> A = B1;

    // Update kernel size for max sigma
    NS = 2 * int(4 * sigma_max + 0.5);

    // Calculate gradients
    Matrix2D<double> B = A;
    Matrix2D<double> GxW1 = ImageProcessor::convMat(B, ImageProcessor::getXGradCore(NS, NS, sigma1));
    Matrix2D<double> GyW1 = ImageProcessor::convMat(B, ImageProcessor::getYGradCore(NS, NS, sigma1));

    // Profile parameters
    int n_sigma = 10;
    int n_myu = 200;
    double sigma_myu = min(sigma1, sigma2);
    int NN = 100;
    int otstup = 15;

    RefinementParameters params;
    params.A = A;
    params.B01 = B01;
    params.NN = NN;
    params.NX = NX;
    params.NY = NY;
    params.n_myu = n_myu;
    params.n_sigma = n_sigma;
    params.otstup = otstup;
    params.sigma0 = sigma0;
    params.sigma01 = sigma01;
    params.sigma1 = sigma1;
    params.sigma2 = sigma2;

    // Thread-safe storage for results
    QVector<QPointF> refinedPositions;
    refinedPositions.resize(selectedEdgePoints.size());

    QVector<double> shiftsX, shiftsY;
    shiftsX.resize(selectedEdgePoints.size());
    shiftsY.resize(selectedEdgePoints.size());

    for (int pointIdx = 0; pointIdx < selectedEdgePoints.size(); pointIdx++) {
        int n0 = selectedEdgePoints[pointIdx].x();
        int m0 = selectedEdgePoints[pointIdx].y();

        double ex = GxW1[n0][m0];
        double ey = GyW1[n0][m0];
        double gradMag = sqrt(ex * ex + ey * ey);

        if (gradMag > 1e-6) {
            ex /= gradMag;
            ey /= gradMag;
        } else {
            ex = 1.0;
            ey = 0.0;
        }

        // Create local copy of params for each thread (thread-safe)
        RefinementParameters localParams = params;
        localParams.ex = ex;
        localParams.ey = ey;

        RefinementResult result = ImageProcessor::refineSinglePoint001(n0, m0, localParams);
        double n_new = result.refinedPosition.x();
        double m_new = result.refinedPosition.y();

        // Thread-safe assignment to specific index
        refinedPositions[pointIdx] = QPointF(n_new, m_new);
        shiftsX[pointIdx] = n_new - n0;
        shiftsY[pointIdx] = m_new - m0;

        if ((pointIdx + 1) % 10 == 0) {
            Logger::instance().log(LogLevel::DEBUG, QString("Refined %1/%2 points").arg(pointIdx + 1).arg(selectedEdgePoints.size()));
        }
    }

    // Draw true edge (green)
    for(int i = 0; i < trueEdge.size(); i++) {
        int x = trueEdge[i].x();
        int y = trueEdge[i].y();
        if (x >= 0 && x < trueEdgeImage.width() && y >= 0 && y < trueEdgeImage.height()) {
            trueEdgeImage.setPixel(x, y, qRgb(0, 255, 0));
        }
    }

    // Draw gradient edges (blue)
    for(int i = 0; i < selectedEdgePoints.size(); i++) {
        int x = selectedEdgePoints[i].x();
        int y = selectedEdgePoints[i].y();
        if (x >= 0 && x < gradientEdgeImage.width() && y >= 0 && y < gradientEdgeImage.height()) {
            gradientEdgeImage.setPixel(x, y, qRgb(0, 0, 255));
        }
    }

    // Draw PCR refined edges (red)
    for(int i = 0; i < refinedPositions.size(); i++) {
        int x = static_cast<int>(refinedPositions[i].x() + 0.5);
        int y = static_cast<int>(refinedPositions[i].y() + 0.5);
        if (x >= 0 && x < pcrEdgeImage.width() && y >= 0 && y < pcrEdgeImage.height()) {
            pcrEdgeImage.setPixel(x, y, qRgb(255, 0, 0));
        }
    }

    // Calculate statistics
    double sumShiftX = 0.0, sumShiftY = 0.0;
    double maxShiftX = 0.0, maxShiftY = 0.0;
    double minShiftX = 1e10, minShiftY = 1e10;

    for (int i = 0; i < selectedEdgePoints.size(); i++) {
        double absShiftX = std::abs(shiftsX[i]);
        double absShiftY = std::abs(shiftsY[i]);

        sumShiftX += absShiftX;
        sumShiftY += absShiftY;

        if (absShiftX > maxShiftX) maxShiftX = absShiftX;
        if (absShiftY > maxShiftY) maxShiftY = absShiftY;

        if (absShiftX < minShiftX) minShiftX = absShiftX;
        if (absShiftY < minShiftY) minShiftY = absShiftY;
    }

    double avgShiftX = (selectedEdgePoints.size() > 0) ? sumShiftX / selectedEdgePoints.size() : 0;
    double avgShiftY = (selectedEdgePoints.size() > 0) ? sumShiftY / selectedEdgePoints.size() : 0;

    // Show statistics
    Logger::instance().log(LogLevel::INFO, "=== REFINEMENT SUMMARY ===");
    Logger::instance().log(LogLevel::INFO, QString("Total points processed: %1").arg(refinedPositions.size()));
    Logger::instance().log(LogLevel::INFO, QString("Average absolute shift: ΔX = %1 pixels, ΔY = %2 pixels").arg(avgShiftX, 0, 'f', 3).arg(avgShiftY, 0, 'f', 3));
    Logger::instance().log(LogLevel::INFO, QString("Maximum absolute shift: ΔX = %1 pixels, ΔY = %2 pixels").arg(maxShiftX, 0, 'f', 3).arg(maxShiftY, 0, 'f', 3));
    Logger::instance().log(LogLevel::INFO, QString("Minimum absolute shift: ΔX = %1 pixels, ΔY = %2 pixels").arg(minShiftX, 0, 'f', 3).arg(minShiftY, 0, 'f', 3));

    // Calculate average shift magnitude
    double avgMagnitude = 0.0;
    double maxMagnitude = 0.0;
    double minMagnitude = 1e10;
    for (int i = 0; i < selectedEdgePoints.size(); i++) {
        double magnitude = std::sqrt(shiftsX[i] * shiftsX[i] + shiftsY[i] * shiftsY[i]);
        avgMagnitude += magnitude;
        if (magnitude > maxMagnitude) maxMagnitude = magnitude;
        if (magnitude < minMagnitude) minMagnitude = magnitude;
    }
    avgMagnitude /= selectedEdgePoints.size();

    Logger::instance().log(LogLevel::INFO, QString("Average shift magnitude: %1 pixels").arg(avgMagnitude, 0, 'f', 3));
    Logger::instance().log(LogLevel::INFO, QString("Maximum shift magnitude: %1 pixels").arg(maxMagnitude, 0, 'f', 3));
    Logger::instance().log(LogLevel::INFO, QString("Minimum shift magnitude: %1 pixels").arg(minMagnitude, 0, 'f', 3));

    // Show final images
    ImageShowcaseWidget* trueEdgeWidget = new ImageShowcaseWidget();
    trueEdgeWidget->setImage(trueEdgeImage);
    trueEdgeWidget->setWindowTitle("Истинная граница (зелёный)");
    trueEdgeWidget->resize(400, 400);
    trueEdgeWidget->show();

    ImageShowcaseWidget* gradientEdgeWidget = new ImageShowcaseWidget();
    gradientEdgeWidget->setImage(gradientEdgeImage);
    gradientEdgeWidget->setWindowTitle("Градиентная граница (синий)");
    gradientEdgeWidget->resize(400, 400);
    gradientEdgeWidget->show();

    ImageShowcaseWidget* pcrEdgeWidget = new ImageShowcaseWidget();
    pcrEdgeWidget->setImage(pcrEdgeImage);
    pcrEdgeWidget->setWindowTitle("Граница ПЦР (красный)");
    pcrEdgeWidget->resize(400, 400);
    pcrEdgeWidget->show();

    // Combined image with all three
    QImage combinedImage = m_dataModel->originalImage();

    for(int i = 0; i < trueEdge.size(); i++) {
        int x = trueEdge[i].x();
        int y = trueEdge[i].y();
        if (x >= 0 && x < combinedImage.width() && y >= 0 && y < combinedImage.height()) {
            combinedImage.setPixel(x, y, qRgb(0, 255, 0));
        }
    }

    for(int i = 0; i < selectedEdgePoints.size(); i++) {
        int x = selectedEdgePoints[i].x();
        int y = selectedEdgePoints[i].y();
        if (x >= 0 && x < combinedImage.width() && y >= 0 && y < combinedImage.height()) {
            // Don't overwrite green
            if (qGreen(combinedImage.pixel(x, y)) != 255) {
                combinedImage.setPixel(x, y, qRgb(0, 0, 255));
            }
        }
    }

    for(int i = 0; i < refinedPositions.size(); i++) {
        int x = static_cast<int>(refinedPositions[i].x() + 0.5);
        int y = static_cast<int>(refinedPositions[i].y() + 0.5);
        if (x >= 0 && x < combinedImage.width() && y >= 0 && y < combinedImage.height()) {
            combinedImage.setPixel(x, y, qRgb(255, 0, 0));
        }
    }

    ImageShowcaseWidget* combinedWidget = new ImageShowcaseWidget();
    combinedWidget->setImage(combinedImage);
    combinedWidget->setWindowTitle("Сравнение: зелёный — истина, синий — градиент, красный — ПЦР");
    combinedWidget->resize(400, 400);
    combinedWidget->show();

    Logger::instance().stopTimer("Edge Refinement Test 001");
    Logger::instance().log(LogLevel::INFO, "=== Edge Refinement Test Completed ===");
}


QImage MainWindow::drawProfileLineOnImage(const QImage& image, const QPoint& p1, const QPoint& p2, const QVector<double>& profile)
{
    QImage result = image.copy();

    // Если изображение цветное, конвертируем в оттенки серого для лучшей видимости
    if (result.format() != QImage::Format_Grayscale8) {
        // Создаем копию и рисуем поверх
    }

    QPainter painter(&result);
    painter.setPen(QPen(Qt::red, 1, Qt::SolidLine));

    // Рисуем линию между точками
    painter.drawLine(p1, p2);
    painter.end();

    return result;
}