#include "dialog.h"
#include "ui_dialog.h"

Dialog::Dialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Dialog)
{
    ui->setupUi(this);
    value = 0.;
}

Dialog::~Dialog()
{
    delete ui;
}

void Dialog::setPlaceholderText(QString newText)
{
    ui->lineEdit->setPlaceholderText(newText);
}

void Dialog::on_lineEdit_editingFinished()
{
    value = ui->lineEdit->text().toDouble();
}

