// dialog.cpp
#include "dialog.h"
#include "ui_dialog.h"

Dialog::Dialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Dialog)
{
    ui->setupUi(this);

    // Set up validator for double input / Ustanavlivaem validator dlya vvoda double
    m_validator = new QDoubleValidator(this);
    m_validator->setDecimals(6);
    m_validator->setRange(-1e6, 1e6);
    m_validator->setLocale(QLocale::c());
    ui->lineEdit->setValidator(m_validator);

    // Connect button box signals / Podklyuchayem signaly knopok
    connect(ui->buttonBox, &QDialogButtonBox::accepted, this, &Dialog::on_buttonBox_accepted);
    connect(ui->buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);
}

Dialog::~Dialog()
{
    delete ui;
}

void Dialog::setPlaceholderText(const QString& newText)
{
    ui->lineEdit->setPlaceholderText(newText);
}

void Dialog::setDefaultValue(double defaultValue)
{
    m_value = defaultValue;
    ui->lineEdit->setText(QString::number(defaultValue));
}

void Dialog::on_lineEdit_editingFinished()
{
    bool ok;
    double val = ui->lineEdit->text().toDouble(&ok);
    if (ok) {
        m_value = val;
    }
}

void Dialog::on_buttonBox_accepted()
{
    // Ensure value is captured even if editing didn't finish
    // Obezpechivayem, chto znacheniye zakhvacheno, dazhe yesli redaktirovaniye ne zavershilos'
    on_lineEdit_editingFinished();
    accept();
    emit setSigma();
}
