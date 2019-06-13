#ifndef RADARCHARACTERISTICS_H
#define RADARCHARACTERISTICS_H

#include <QVector>
#include <cmath>

// Класс для расчета характеристик РЛМ

class RadarCharacteristics
{
private:
    // Общие константы
    const double RAD = M_PI/180; // Константа для перевода из радиан в градусы
    const double C = 299792458; // Скорость света
public:
    // ДН по азимуту
    QVector<double> DN_beta(double Nx, double dx, double f, double beta_target, const QVector<double> &beta);
    // ДН по углу места
    QVector<double> DN_eps(double Ny, double dy, double f, double eps_a,
                           double eps_target, double Hgc, bool reflection, const QVector<double> &eps);

    // Пеленгационная характеристика
    QVector<double> Peleng(double Ny, double dy, double f, double eps_a,
                           double Hgc, bool reflection, const QVector<double> &eps);

    // Проводка по высоте
    // Дальности
    QVector<double> Distances(double H0, const QVector<double> &eps);
    // Высоты
    QVector<double> Heights(const QVector<double> &eps, const QVector<double> &eps_real, const QVector<double> &distance);
};

#endif // RADARCHARACTERISTICS_H
