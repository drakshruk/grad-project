// dialog.h
#ifndef DIALOG_H
#define DIALOG_H

#include <QDialog>
#include <QDoubleValidator>

namespace Ui {
class Dialog;
}

/**
 * EN: Simple input dialog for numeric values (sigma, radius, etc.)
 * RU: Prostoy dialog vvoda dlya chislovykh znacheniy (sigma, radius, i t.d.)
 */
class Dialog : public QDialog
{
    Q_OBJECT

public:
    explicit Dialog(QWidget *parent = nullptr);
    ~Dialog();

    /**
     * EN: Gets the entered value
     * RU: Poluchayet vvedennoye znacheniye
     */
    double getValue() const { return m_value; }

    /**
     * EN: Sets the placeholder text in the input field
     * RU: Ustanavlivayet tekst-podskazku v pole vvoda
     */
    void setPlaceholderText(const QString& newText);

    /**
     * EN: Sets the default value
     * RU: Ustanavlivayet znacheniye po-umolchaniyu
     */
    void setDefaultValue(double defaultValue);

signals:
    /**
     * EN: Emitted when sigma value is confirmed
     * RU: Ispuskayetsya pri podtverzhdenii znacheniya sigma
     */
    void setSigma();

private slots:
    void on_lineEdit_editingFinished();
    void on_buttonBox_accepted();

private:
    Ui::Dialog *ui;
    double m_value = 0.0;
    QDoubleValidator* m_validator = nullptr;
};

#endif // DIALOG_H