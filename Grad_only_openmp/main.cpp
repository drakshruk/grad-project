#include "mainwindow.h"
#include <QApplication>
#include <omp.h>

int main(int argc, char *argv[])
{
    omp_set_num_threads(12);  // Фиксированное значение

    // int numThreads = omp_get_max_threads();
    // omp_set_num_threads(numThreads);

    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}