# Grad — Image Processing & Edge Refinement Toolkit

[![Qt](https://img.shields.io/badge/Qt-5.15%2B%20%7C%206.x-41CD52?logo=qt&logoColor=white)](https://www.qt.io/)
[![C++](https://img.shields.io/badge/C%2B%2B-17-00599C?logo=c%2B%2B&logoColor=white)](https://en.cppreference.com/w/cpp/17)
[![OpenMP](https://img.shields.io/badge/OpenMP-enabled-EE4C2C?logo=openmp&logoColor=white)](https://www.openmp.org/)
[![QCustomPlot](https://img.shields.io/badge/QCustomPlot-2.x-1E88E5)](https://www.qcustomplot.com/)

A Qt-based desktop application for **image processing**, **gradient/Laplacian edge detection**, and **sub-pixel edge refinement** using a custom algorithm.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Tech Stack](#tech-stack)
- [Project Structure](#project-structure)
- [Development History](#development-history)
- [Building & Running](#building--running)
- [Usage](#usage)
- [Algorithm Notes](#algorithm-notes)
- [License](#license)
- [Author](#author)

---

## Overview

**Grad** is a research/engineering GUI tool built around a custom image-processing pipeline. It started as a simple gradient-based edge detector and gradually grew into a full playground for testing convolution kernels, edge-detection operators, radiography-style simulations, and a proprietary sub-pixel edge-refinement algorithm.

The application is written in **C++17** using **Qt Widgets**, accelerated with **OpenMP**, and visualizes profiles and signals with **QCustomPlot**.

> **Note:** The project is primarily an experimental/research tool — internal APIs and file formats may change.

---

## Features

### Image I/O
- Open / save images (PNG, JPG, BMP, TIFF, XPM).
- Built-in sample generators: `two hollows`, `two hollows big`, synthetic radiographic cylinder with cavity.

### Filtering & Edge Detection
- **Gaussian blur** with configurable σ and kernel radius.
- **Sharpen** (3×3 Laplacian-based kernel).
- **Gaussian edge detection** (Difference of Gaussians + zero-crossing).
- **Laplacian edge detection** (LoG kernel + zero-crossing).
- **Gradient X / Y** via Gaussian derivatives; saves `XGradImage.png`, `YGradImage.png`.

### Profiles
- Vertical profile through the image center.
- **Profile between two points** — click two points, get the intensity profile along the line (Bresenham or bilinear interpolation).
- Multi-profile comparison with labels and axis annotations.
- **Zero-crossing visualization** on Laplacian profiles.

### Sub-pixel Edge Refinement (PCR)
- Interactive **edge selection** (flood-fill of connected edge pixels).
- Per-point **Profile Correlation Refinement** — searches for the optimal shift and scale of a blurred profile relative to a reference profile.
- Refined points reported with sub-pixel coordinates.
- Statistics: mean / max / min shift, shift magnitude, residual.
- Side-by-side comparison: *true edge* / *gradient edge* / *PCR-refined edge*.

### Image Calculator
- Element-wise operations between two images:
  `Add`, `Subtract`, `Multiply`, `Divide`, `AND`, `OR`, `XOR`, `Min`, `Max`, `Average`, `Difference`, `Copy`, `Transparent-zero`.
- Result can go to the main window or a new one.

### Statistics & Visualization
- Mean, median, min, max, range, variance, std-dev, skewness, kurtosis, entropy.
- 256-bin histogram.
- Export statistics to a text file.

### Logging
- Singleton `Logger` writing to `./logs/app_log_<timestamp>.txt`.
- Colorized in-app log window (`LoggerWidget`).
- Per-method timers via `Logger::startTimer` / `stopTimer` and an RAII `ScopedTimer`.
- On shutdown, the session log is copied to `./logs/last_session.log`.

### Benchmarking
- `benchmarkConvolution()` and `benchmarkRefinement()` — compare 1/2/4/8 OpenMP threads and report speedup / efficiency.

---

## Tech Stack

| Component | Purpose |
|---|---|
| **C++17** | Language |
| **Qt 5.15+ / 6.x** (Widgets, PrintSupport) | GUI |
| **OpenMP** | Parallel convolution / edge detection / refinement |
| **QCustomPlot 2.x** | Profile & signal plotting |
| **qmake** (`.pro`) | Build system |
| **MinGW / MSVC / GCC / Clang** | Compilers |

---

## Project Structure

```
Grad/
├── Grad.pro                     # qmake project file
├── main.cpp                     # Entry point (sets OpenMP threads)
├── mainwindow.{h,cpp,ui}        # Main window, actions, feature wiring
├── appdatamodel.{h,cpp}         # Central data model (images, matrices, params)
├── imageprocessingtypes.h       # Matrix aliases, structs (Refinement*, Profile*, App_Stats)
├── imageprocessor.{h,cpp}       # Core algorithms: convolution, kernels, edges, refinement
├── matrixlambdas.h              # Functors for element-wise matrix operations
├── imagecalculator.{h,cpp,ui}   # Image arithmetic UI
├── imageshowcasewidget.{h,cpp}  # Clickable image display widget
├── graphwidget.{h,cpp,ui}       # QCustomPlot wrapper (profiles, zero-crossings)
├── dialog.{h,cpp,ui}            # Numeric input dialog
├── loggerwidget.{h,cpp,ui}      # Singleton logger + log window
├── qcustomplot.{h,cpp}          # Vendored QCustomPlot
└── logs/                        # Runtime log output (auto-created)
```

### Key modules

| Module | Role |
|---|---|
| `ImageProcessor` | Stateless namespace with all algorithms: `convMat`, `getGauss`, `getXGradCore`, `getYGradCore`, `getLapl`, `findEdges`, `refineSinglePoint001/002`, `buildProfileBetweenPoints`, `computeImageStatistics`, … |
| `AppDataModel` | Single source of truth for current/original image, image matrix, attribute matrix, selected edge, parameters (σ, radius, …). Emits Qt signals on change. |
| `Logger` / `LoggerWidget` | Thread-safe singleton with buffered logs, file logging, timers, and a GUI window. |
| `GraphWidget` | Wraps `QCustomPlot` and adds helpers: `plotMultipleProfiles`, `plotMultipleProfilesWithZeroCrossings`. |
| `ImageShowcaseWidget` | Renders a `QPixmap` with aspect-ratio preservation and converts widget ↔ image coordinates on click. |

---

## Development History

The project evolved incrementally over ~3+ months:

1. **Gradient-based edge detection** — the original core. Gaussian / Gaussian-derivative convolutions, LoG, DoG, zero-crossing detection.
2. **Edge selection UI** — flood-fill selection of connected edge pixels from a click, visualized in cyan.
3. **Profile visualization** — QCustomPlot integration, single & multi-profile plots, profile between two user-clicked points.
4. **Profile Correlation Refinement (PCR)** — the main algorithmic contribution: compare a reference (sharp) Laplacian profile against a blurred Laplacian profile, find the shift & scale that minimizes the normalized residual, and convert that shift back to sub-pixel edge displacement.
5. **Supporting infrastructure** — image calculator, statistics, sample generators, benchmarking harness, logging subsystem, and a full QWidget UI.

Mathematics (kernels, normalizations, residual formulas) was added *ad hoc* as the algorithms were prototyped — so some naming and comments are intentionally verbose (bilingual EN/RU) to keep the code self-documenting.

---

## Building & Running
 (not tested)
### Requirements
- Qt **5.15** or **6.x** (`core`, `gui`, `widgets`, `printsupport`)
- A C++17 compiler (MinGW, MSVC 2019+, GCC 9+, Clang 10+)
- OpenMP runtime (bundled with GCC/MinGW; `/openmp` on MSVC)
- qmake (comes with Qt)

### Build (qmake)

```bash
git clone https://github.com/<your-username>/Grad.git
cd Grad
mkdir build && cd build
qmake ../Grad.pro
make -j$(nproc)         # Linux / macOS
# or
mingw32-make -j8        # MinGW on Windows
# or open Grad.pro in Qt Creator and press Run
```

### OpenMP notes
The `.pro` file already enables OpenMP for:
- `win32-g++` → `-fopenmp`
- `win32-msvc*` → `/openmp`
- `unix:!macx` → `-fopenmp`

`main.cpp` currently hard-codes `omp_set_num_threads(12)` — change this to match your CPU:

```cpp
omp_set_num_threads(12);  // ← adjust
// or, to use the hardware maximum:
// omp_set_num_threads(omp_get_max_threads());
```

### Release build
The project ships with optimization flags:

```
QMAKE_CXXFLAGS_RELEASE += -O3 -march=native
QMAKE_CXXFLAGS_DEBUG   += -O0 -g
```

`-march=native` is dropped automatically if you switch to MSVC.

---

## Usage

1. **Launch** `Grad`.
2. **File → open file** or **Samples → two hollows** to load an image.
3. **Math menu:**
   - *gaussian blur* — enter σ, then kernel radius.
   - *gaussian edge detection* / *Laplacian edge detection* — enter σ; the result is drawn and the app enters **edge-selection mode**.
   - *gradient X and Y* — writes `XGradImage.png`, `YGradImage.png`.
   - *sharpen* — 3×3 unsharp kernel.
4. **Click on a detected edge** in the image window to select it (flood-fill of connected pixels).
5. **Tools → test_001** or **test_002** to run PCR refinement on the selected edge / on a synthetic two-hollow sample. Results are shown in separate windows (true / gradient / PCR).
6. **Tools → draw profile** — vertical profile through the center.
7. **Tools → build profile** — click two points to plot the profile along the line.
8. **Tools → image calculator** — pick two images and an operation.
9. **Tools → show statistics / export statistics** — analyze the current image.
10. **Tools → show log window** — inspect the runtime log.

---

## Algorithm Notes

### Profile Correlation Refinement (PCR) — in short

Given an initial edge point `(n₀, m₀)` and the gradient direction `(ex, ey)`:

1. Sample the LoG response of two images (reference `A`, blurred `B01`) along the gradient direction, producing profiles `prof1(μ)` and `prof2(μ)`.
2. Find the zero-crossings of each profile (`nL_Zero`, `nR_Zero`).
3. Detect the interval `[XL1, XR1]` around the zero-crossing by scanning for plateau (5 consecutive increases/decreases).
4. Reinterpolate `prof1` over `[XL1, XR1]` to a fixed length `NN`.
5. Search for the best `[XL2, XR2]` in the neighborhood (`±otstup`) of `[XL1, XR1]` that minimizes the **normalized RMS residual** between the two profiles.
6. Convert the optimal bounds into a **shift** `FFF1_1` and a **scale** `FFF1_0`, and compute the refined point as
   `(n_new, m_new) = (n₀ + ex · FFF1_1, m₀ + ey · FFF1_1)`.

The second variant (`refineSinglePoint002`) replaces the exhaustive `(2·otstup+1)²` search with a **gradient-descent search** and adds an `exp()` LUT, giving a noticeable speedup for large `otstup`.

### Kernels
- `getGauss(iX, iY, σ)` — separable Gaussian (not normalized; the `convMat` divides by the kernel sum).
- `getXGradCore`, `getYGradCore` — Gaussian derivatives.
- `getLapl(iX, iY, σ)` — Laplacian of Gaussian: `(r²/σ² − 2)·exp(−r²/2σ²)`.

### Parallelism
The heavy loops (`convMat`, `findEdges`, `refineSinglePoint001/002`) are annotated with `#pragma omp parallel for` (some are currently commented out — uncomment to enable).

---

## License

This project is provided **as-is** for research and educational purposes.
QCustomPlot is licensed under the **GPLv3** — see `qcustomplot.h` for details.
If you plan to redistribute, make sure your Qt license and QCustomPlot license are compatible with your intended use.

---

## Author

**Yaroslav Konovalov**
GitHub: [@drakshruk](https://github.com/drakshruk)

---

<details>
<summary><b>🇷🇺 Русская версия</b></summary>

# Grad — инструмент обработки изображений и уточнения границ

Qt-приложение для **обработки изображений**, **градиентной и лапласианской детекции границ** и **субпиксельного уточнения положения границ** с помощью алгоритма ПЦР (поиска центров растяжения).

---

## Содержание

- [Обзор](#обзор)
- [Возможности](#возможности)
- [Стек технологий](#стек-технологий)
- [Структура проекта](#структура-проекта)
- [История разработки](#история-разработки)
- [Сборка и запуск](#сборка-и-запуск)
- [Использование](#использование)
- [Заметки об алгоритмах](#заметки-об-алгоритмах)
- [Лицензия](#лицензия)
- [Автор](#автор)

---

## Обзор

**Grad** — это GUI-инструмент, выросший из простого градиентного детектора границ в полноценную площадку для экспериментов со свёртками, операторами выделения границ, симуляцией радиографических объектов и собственным алгоритмом субпиксельного уточнения.

Проект написан на **C++17** с использованием **Qt Widgets**, ускорен **OpenMP**, а визуализация профилей и сигналов выполнена через **QCustomPlot**.

> **Внимание:** проект носит исследовательский характер — внутренние API и форматы могут меняться.

---

## Возможности

### Ввод/вывод изображений
- Открытие / сохранение PNG, JPG, BMP, TIFF, XPM.
- Встроенные генераторы примеров: `two hollows`, `two hollows big`, синтетический радиографический цилиндр с полостью.

### Фильтрация и детекция границ
- **Гауссово размытие** с настраиваемыми σ и радиусом ядра.
- **Sharpen** — 3×3 ядро на основе лапласиана.
- **Gaussian edge detection** — Difference of Gaussians + zero-crossing.
- **Laplacian edge detection** — ядро LoG + zero-crossing.
- **Gradient X / Y** через производные Гаусса; сохраняет `XGradImage.png`, `YGradImage.png`.

### Профили
- Вертикальный профиль через центр изображения.
- **Профиль между двумя точками** — кликните две точки и получите профиль вдоль линии (Брезенхэм или билинейная интерполяция).
- Мультипрофильное сравнение с подписями и осями.
- **Визуализация zero-crossing** на профилях Лапласиана.

### Субпиксельное уточнение границ (PCR)
- Интерактивный **выбор границы** (flood-fill связных краевых пикселей).
- **Profile Correlation Refinement** для каждой точки — поиск оптимального сдвига и масштаба размытого профиля относительно опорного.
- Уточнённые точки с субпиксельными координатами.
- Статистика: средний / макс / мин сдвиг, модуль сдвига, невязка.
- Сравнение: *истинная граница* / *градиентная* / *PCR*.

### Калькулятор изображений
- Поэлементные операции между двумя изображениями:
  `Add`, `Subtract`, `Multiply`, `Divide`, `AND`, `OR`, `XOR`, `Min`, `Max`, `Average`, `Difference`, `Copy`, `Transparent-zero`.
- Результат — в главное окно или в новое.

### Статистика и визуализация
- Среднее, медиана, min, max, размах, дисперсия, СКО, асимметрия, эксцесс, энтропия.
- Гистограмма на 256 бинов.
- Экспорт статистики в текстовый файл.

### Логирование
- Singleton `Logger` пишет в `./logs/app_log_<timestamp>.txt`.
- Цветное окно лога в приложении (`LoggerWidget`).
- Таймеры методов через `Logger::startTimer` / `stopTimer` и RAII-обёртку `ScopedTimer`.
- При выходе сессионный лог копируется в `./logs/last_session.log`.

### Бенчмарки
- `benchmarkConvolution()` и `benchmarkRefinement()` — сравнение 1/2/4/8 потоков OpenMP, вывод ускорения и эффективности.

---

## Стек технологий

| Компонент | Назначение |
|---|---|
| **C++17** | Язык |
| **Qt 5.15+ / 6.x** (Widgets, PrintSupport) | GUI |
| **OpenMP** | Параллельные свёртки / детекция / уточнение |
| **QCustomPlot 2.x** | Графики профилей и сигналов |
| **qmake** (`.pro`) | Система сборки |
| **MinGW / MSVC / GCC / Clang** | Компиляторы |

---

## Структура проекта

```
Grad/
├── Grad.pro                     # qmake-проект
├── main.cpp                     # Точка входа (устанавливает число потоков OpenMP)
├── mainwindow.{h,cpp,ui}        # Главное окно, экшены, связывание функций
├── appdatamodel.{h,cpp}         # Центральная модель данных (изображения, матрицы, параметры)
├── imageprocessingtypes.h       # Алиасы матриц, структуры (Refinement*, Profile*, App_Stats)
├── imageprocessor.{h,cpp}       # Ядро алгоритмов: свёртки, ядра, границы, уточнение
├── matrixlambdas.h              # Функторы для поэлементных матричных операций
├── imagecalculator.{h,cpp,ui}   # UI калькулятора изображений
├── imageshowcasewidget.{h,cpp}  # Виджет отображения изображения с обработкой кликов
├── graphwidget.{h,cpp,ui}       # Обёртка над QCustomPlot (профили, zero-crossing)
├── dialog.{h,cpp,ui}            # Диалог ввода числа
├── loggerwidget.{h,cpp,ui}      # Singleton-логгер + окно лога
├── qcustomplot.{h,cpp}          # Вендоренный QCustomPlot
└── logs/                        # Логи во время работы (создаётся автоматически)
```

### Ключевые модули

| Модуль | Роль |
|---|---|
| `ImageProcessor` | Namespace без состояния: `convMat`, `getGauss`, `getXGradCore`, `getYGradCore`, `getLapl`, `findEdges`, `refineSinglePoint001/002`, `buildProfileBetweenPoints`, `computeImageStatistics`, … |
| `AppDataModel` | Единый источник истины: текущее/оригинальное изображение, матрицы, выбранная граница, параметры (σ, радиус, …). Шлёт Qt-сигналы при изменениях. |
| `Logger` / `LoggerWidget` | Потокобезопасный singleton с буфером логов, записью в файл, таймерами и GUI-окном. |
| `GraphWidget` | Оборачивает `QCustomPlot` и добавляет `plotMultipleProfiles`, `plotMultipleProfilesWithZeroCrossings`. |
| `ImageShowcaseWidget` | Рисует `QPixmap` с сохранением пропорций и конвертирует координаты виджет ↔ изображение по клику. |

---

## История разработки

Проект развивался итеративно в течение ~3+ месяцев:

1. **Градиентная детекция границ** — исходное ядро: свёртки Гаусса и его производных, LoG, DoG, zero-crossing.
2. **UI выбора границы** — flood-fill связных краевых пикселей от клика, подсветка голубым.
3. **Визуализация профилей** — интеграция QCustomPlot, одиночные и мультипрофили, профиль между двумя кликами.
4. **Profile Correlation Refinement (PCR)** — ключевой алгоритм: сравнение опорного (резкого) и размытого лапласианских профилей, поиск сдвига и масштаба, минимизирующих нормированную невязку, и перевод сдвига обратно в субпиксельное смещение границы.
5. **Инфраструктура** — калькулятор изображений, статистика, генераторы примеров, бенчмарки, подсистема логирования, полноценный QWidget-интерфейс.

Математика (ядра, нормировки, формулы невязок) добавлялась *по мере необходимости* — поэтому комментарии намеренно подробные и двуязычные (EN/RU), чтобы код оставался самодокументируемым.

---

## Сборка и запуск
(не протестировано)
### Требования
- Qt **5.15** или **6.x** (`core`, `gui`, `widgets`, `printsupport`)
- Компилятор с C++17 (MinGW, MSVC 2019+, GCC 9+, Clang 10+)
- Runtime OpenMP (идёт с GCC/MinGW; `/openmp` для MSVC)
- qmake (идёт с Qt)

### Сборка (qmake)

```bash
git clone https://github.com/<your-username>/Grad.git
cd Grad
mkdir build && cd build
qmake ../Grad.pro
make -j$(nproc)         # Linux / macOS
# или
mingw32-make -j8        # MinGW на Windows
# или откройте Grad.pro в Qt Creator и нажмите Run
```

### OpenMP
В `.pro` уже включён OpenMP для:
- `win32-g++` → `-fopenmp`
- `win32-msvc*` → `/openmp`
- `unix:!macx` → `-fopenmp`

`main.cpp` жёстко задаёт `omp_set_num_threads(12)` — поменяйте под свой CPU:

```cpp
omp_set_num_threads(12);  // ← измените
// или, чтобы использовать максимум:
// omp_set_num_threads(omp_get_max_threads());
```

### Release-сборка
В проекте уже есть флаги оптимизации:

```
QMAKE_CXXFLAGS_RELEASE += -O3 -march=native
QMAKE_CXXFLAGS_DEBUG   += -O0 -g
```

`-march=native` автоматически отключается при сборке MSVC.

---

## Использование

1. **Запустите** `Grad`.
2. **File → open file** или **Samples → two hollows**.
3. **Меню Math:**
   - *gaussian blur* — введите σ, затем радиус ядра.
   - *gaussian edge detection* / *Laplacian edge detection* — введите σ; результат отрисуется, приложение войдёт в **режим выбора границы**.
   - *gradient X and Y* — сохранит `XGradImage.png`, `YGradImage.png`.
   - *sharpen* — ядро 3×3.
4. **Кликните по найденной границе** — flood-fill выберет связные пиксели.
5. **Tools → test_001** / **test_002** — запуск PCR-уточнения для выбранной границы / синтетического примера. Результаты показываются в отдельных окнах (истинная / градиентная / PCR).
6. **Tools → draw profile** — вертикальный профиль через центр.
7. **Tools → build profile** — кликните две точки, чтобы построить профиль вдоль линии.
8. **Tools → image calculator** — выберите два изображения и операцию.
9. **Tools → show statistics / export statistics** — анализ текущего изображения.
10. **Tools → show log window** — окно лога.

---

## Заметки об алгоритмах

### Profile Correlation Refinement (PCR) — кратко

Дано: начальная точка границы `(n₀, m₀)` и направление градиента `(ex, ey)`.

1. Сэмплируем LoG-отклик двух изображений (опорного `A` и размытого `B01`) вдоль градиента → профили `prof1(μ)`, `prof2(μ)`.
2. Находим zero-crossing каждого профиля (`nL_Zero`, `nR_Zero`).
3. Определяем интервал `[XL1, XR1]` вокруг zero-crossing сканированием плато (5 подряд возрастаний/убываний).
4. Переинтерполируем `prof1` на `[XL1, XR1]` до фиксированной длины `NN`.
5. Ищем лучший `[XL2, XR2]` в окрестности (`±otstup`) от `[XL1, XR1]`, минимизируя **нормированную RMS-невязку** между профилями.
6. Переводим оптимальные границы в **сдвиг** `FFF1_1` и **масштаб** `FFF1_0` и вычисляем уточнённую точку:
   `(n_new, m_new) = (n₀ + ex · FFF1_1, m₀ + ey · FFF1_1)`.

Вариант `refineSinglePoint002` заменяет полный перебор `(2·otstup+1)²` на **градиентный спуск** и добавляет LUT для `exp()`, что даёт заметное ускорение при больших `otstup`.

### Ядра
- `getGauss(iX, iY, σ)` — разделяемое гауссово (не нормировано; `convMat` делит на сумму ядра).
- `getXGradCore`, `getYGradCore` — производные Гаусса.
- `getLapl(iX, iY, σ)` — LoG: `(r²/σ² − 2)·exp(−r²/2σ²)`.

### Параллелизм
Тяжёлые циклы (`convMat`, `findEdges`, `refineSinglePoint001/002`) помечены `#pragma omp parallel for` (часть закомментирована — раскомментируйте для включения).

---

## Лицензия

Проект предоставляется **как есть** для исследовательских и образовательных целей.
QCustomPlot распространяется под **GPLv3** — см. `qcustomplot.h`.
При распространении убедитесь, что лицензии Qt и QCustomPlot совместимы с вашим сценарием.

---

## Автор

**Ярослав Коновалов**
GitHub: [@drakshruk](https://github.com/drakshruk)

</details>
