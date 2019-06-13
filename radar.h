/****************************************************************************
**           Author: Artem Kudryavtsev                                    **
**             Date: 24.04.19                                             **
**          Version: 1.1.0                                                **
****************************************************************************/

#ifndef RADAR_H
#define RADAR_H

#include <QMainWindow>
#include <cmath>
#include <qcustomplot.h>

namespace Ui {
class Radar;
}

// Главный классс приложения

class Radar : public QMainWindow
{
    Q_OBJECT
private:
    // Параметры для ДН по азимуту
    double Nx; // Количество столбцов АР
    double dx; // Расстояние между элементами
    double f_beta; // Частота
    double beta_target; // Угол цели (на который настроен пространственный фильтр)
    QVector<double> beta; // Массив азимутов
    QVector<double> DN_beta; // Значения амплитуд ДН

    // Параметры для ДН по углу места
    double Ny; // Количество строк АР
    double dy; // Расстояние между элементами
    double f_eps; // Частота
    double eps_a; // Наклон АР
    double Hgc; // Высота подъема геометрического центра АР
    double eps_target; // Угол цели (на который настроен пространственный фильтр)
    QVector<double> eps; // Массив углов места
    QVector<double> DN_eps; // Значения амплитуд ДН

    // Параметры для пеленгационной характеристики
    double Ny_pel; // Количество строк
    double dy_pel; // Расстояние между элементами
    double f_eps_pel; // Частота
    double eps_a_pel; // Наклон АР
    double Hgc_pel; // Высота подъема геометрического центра АР
    QVector<double> eps_orig_pel; // Вектор истинных углов места
    QVector<double> eps_real_pel; // Вектор измеренных углов места

    // Параметры для проводки по высоте
    double Ny_prov; // Количество строк
    double dy_prov; // Расстояние между элементами
    double f_eps_prov; // Частота
    double eps_a_prov; // Наклон АР
    double Hgc_prov; // Высота подъема геометрического центра АР
    double H0_prov; // Высота полета
    QVector<double> eps_orig_prov; // Вектор истинных углов места
    QVector<double> eps_real_prov; // Вектор измеренных углов места
    QVector<double> distance_prov; // Дальности
    QVector<double> height_prov; // Высоты

public:
    explicit Radar(QWidget *parent = 0);
    ~Radar();

private slots:
    // ДН по азимуту
    void on_buttonDefaultBeta_clicked(); // Параметры по умолчанию
    void on_buttonPlotBeta_clicked(); // Построение графика
    void on_buttonSaveBeta_clicked(); // Сохранение в файл

    // ДН по углу места
    void on_buttonDefaultEps_clicked(); // Параметры по умолчанию
    void on_buttonPlotEps_clicked(); // Построение графика
    void on_buttonSaveEps_clicked(); // Сохранение в файл

    // Пеленгационная характеристика
    void on_buttonDefaultPeleng_clicked(); // Параметры по умолчанию
    void on_buttonPlotPeleng_clicked(); // Построение графика
    void on_buttonSavePeleng_clicked(); // Сохранение в файл

    // Проводка по высоте
    void on_buttonDefaultProvodka_clicked(); // Параметры по умолчанию
    void on_buttonPlotProvodka_clicked(); // Построение графика
    void on_buttonSaveProvodka_clicked(); // Сохранение в файл

    // Выход из приложения через меню
    void exit() {
        emit QWidget::close();
    }

private:
    Ui::Radar *ui;
    // Функция построения графика по двум векторам
    void plot(const QVector<double> &x, const QVector<double> &y, QCustomPlot *area);
};

#endif // RADAR_H
