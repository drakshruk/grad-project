#ifndef DIALOG_H
#define DIALOG_H

#include <QDialog>

namespace Ui {
class Dialog;
}

class Dialog : public QDialog
{
    Q_OBJECT

public:
    explicit Dialog(QWidget *parent = nullptr);
    ~Dialog();
    double value;
    void setPlaceholderText(QString newText);
signals:
    void setSigma();

private slots:
    void on_lineEdit_editingFinished();

private:
    Ui::Dialog *ui;
};

#endif // DIALOG_H
