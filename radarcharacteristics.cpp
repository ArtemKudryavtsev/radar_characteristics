#include "radarcharacteristics.h"
#include <cmath>
#include <complex>
#include <QtDebug>
#include <fstream>

// ДН по азимуту
QVector<double> RadarCharacteristics::DN_beta(double Nx, double dx, double f,
                                              double beta_target, const QVector<double> & beta)
{
    // Nx - Количество элементов по x
    // dx - Расстояние между элементами по x
    // f - Частота
    // beta_target - Угол цели (на который настроен пространственный фильтр)
    // beta - Массив азимутов
     double lam = C/f; // Длина волны
     double k = 2 * M_PI / lam; // Волновое число

    QVector<double> DN_beta(beta.size()); // ДН
    int j = 0;
    for (int i = 0; i < beta.size(); i++)
    {
        std::complex<double> DN_i(0.0, 0.0);
        for (int n = 0; n < Nx; n++)
        {
            // Диаграмма направленности элемента
            double DN_elem = cos(beta[i] * RAD);
            // Сигнал прямой
            std::complex<double> S = 1i * k * sin(beta[i] * RAD) * dx * n;
            std::complex<double> S0 = -1i * k * sin(beta_target * RAD) * dx * n;
            // Весовая функция
            double Wn = 1;
            DN_i += DN_elem * Wn * exp(S) * exp(S0);
        }
        DN_beta[j] = abs(DN_i);
        j++;
    }
    return DN_beta;
}

// ДН по углу места
QVector<double> RadarCharacteristics::DN_eps(double Ny, double dy, double f, double eps_a,
                                             double eps_target, double Hgc,
                                             bool reflection, const QVector<double> &eps)
{
    // Ny - Количество элементов по x
    // dy - Расстояние между элементами по x
    // f - Частота
    // eps_a - Наклон АР
    // eps_target - Угол цели (на который настроен пространственный фильтр)
    // Hgc - Высота подъема геоиетрического центра АР
    // reflection - Отражения от земли (есть/нет)
    // eps - Массив углов места

    double lam = C/f; // Длина волны
    double k = 2 * M_PI / lam; // Волновое число

    double sigma = 0.001; // Проводимость
    double dial_pr = 10; // Относительная диэлектрическая проницаемость
    // Относительная комплексная диэлектрическая проницаемость
    std::complex<double> dial_pr_complex = dial_pr - 1i * 60 * lam * sigma;


    QVector<double> DN_eps(eps.size()); // ДН
    int j = 0;
    for (int i = 0; i < eps.size(); i++)
    {
        // Коэффициент отражения от земной поверхности
        std::complex<double> Kotr(0,0); // Ноль
        if (reflection != false) // Иначе, если отражения есть
            Kotr = (dial_pr_complex * sin(eps[i] * RAD) - sqrt(dial_pr_complex - pow(cos(eps[i] * RAD),2))) /
                    (dial_pr_complex * sin(eps[i] * RAD) + sqrt(dial_pr_complex - pow(cos(eps[i] * RAD),2)));
        std::complex<double> DN_i(0.0, 0.0);
        for (int n = 0; n < Ny; n++)
        {
            double DN_elem_pr = cos((eps[i] - eps_a) * RAD); // Диаграмма направленности элемента
            double DN_elem_otr = cos((eps[i] + eps_a) * RAD); // Диаграмма направленности элемента
            // Сигнал прямой
            std::complex<double> Spr = 1i * k * sin(eps[i] * RAD) * dy * n;
            // Сигнал отраженный
            std::complex<double> Sotr = -1i * k * (2 * (Hgc - cos(eps_a * RAD) * dy * (Ny-1) / 2) *
                                                sin(eps[i] * RAD) + sin((eps[i] + eps_a) * RAD) * dy * n);
            // Вес пространственного фильтра
            std::complex<double> Wn = -1i * k * sin(eps_target * RAD) * dy * n;
            DN_i += exp(Wn) * (DN_elem_pr * exp(Spr) + Kotr * DN_elem_otr * exp(Sotr));
        }
        DN_eps[j] = abs(DN_i);
        j++;
    }
    return DN_eps;
}

// Пеленгационная характеристика
QVector<double> RadarCharacteristics::Peleng(double Ny, double dy, double f, double eps_a,
                                             double Hgc, bool reflection, const QVector<double> & eps)
{
    // Ny - Количество элементов по x
    // dy - Расстояние между элементами по x
    // f - Частота
    // eps_a - Наклон АР
    // Hgc - Высота подъема геометрического центра АР
    // reflection - Отражения от земли (есть/нет)
    // eps - Массив углов места

    double lam = C/f; // Длина волны
    double k = 2 * M_PI / lam; // Волновое число

    double sigma = 0.001; // Проводимость
    double diel = 10; // Относительная диэлектрическая проницаемость
    std::complex<double> diel_compl = diel - 1i * 60 * lam * sigma; // Относительная комплексная диэлектрическая проницаемость

    double eps_start = eps[0]; // Начальный угол сканирования
    double eps_step = eps[1] - eps[0]; // Шаг сканирования
    double eps_stop = eps[eps.size() - 1]; // Конечный угол сканирования
    int eps_size = (eps_stop - eps_start) / eps_step + 1; // Размер массива

    QVector<double> eps_real(eps_size); // Массив измеренных углов места
    for (int i = 0; i < eps_size; i++)
    {
        QVector<double> DN_epsL(eps_size);
        // Коэффициент отражения от земной поверхности (вертикальная поляризация)
        std::complex<double> Kotr(0.0); // Отражений нет
        if (reflection == true) // Отражения есть
        {
            Kotr = (diel_compl * sin(eps[i] * RAD) - sqrt(diel_compl - pow(cos(eps[i] * RAD), 2))) /
                    (diel_compl * sin(eps[i] * RAD) + sqrt(diel_compl - pow(cos(eps[i] * RAD), 2)));
        }
        for (int j = 0; j < eps_size; j++) // Набор пространственных фильтров
        {
            std::complex<double> Y(0.0);
            for (int n = 0 ; n < Ny; n++)
            {
                double DN_elem_eps_pr = pow(cos((eps[i] - eps_a) * RAD), 2); // Диаграмма направленности элемента
                double DN_elem_eps_otr = pow(cos((eps[i] + eps_a) * RAD), 2); // Диаграмма направленности элемента мнимой АР

                std::complex<double> W_n_arg = -1i * k * sin((eps[j] - eps_a) * RAD) * dy * n;
                std::complex<double> W_n = exp(W_n_arg); // Коэффициенты пространственного фильтра

                std::complex<double> Spr_arg = 1i * k * sin((eps[i] - eps_a) * RAD) * dy * n;
                std::complex<double> Spr = exp(Spr_arg); // Прямое колебание

                std::complex<double> Sotr_arg = -1i * k * (2 * (Hgc - cos(eps_a * RAD) * dy * (Ny - 1) / 2) *
                                                           sin(eps[i] * RAD) + sin((eps[i] + eps_a) * RAD) * dy * n);
                std::complex<double> Sotr = exp(Sotr_arg); // Отражённое колебание

                Y = Y + W_n * (Spr * DN_elem_eps_pr  + Sotr * DN_elem_eps_otr * Kotr);
            }
            DN_epsL[j] = abs(Y);
        }

        // Поиск номера элемента с максимальной амплитудой
        // для определения измеренного угла места
        int max_index = 0;
        double max_DN = DN_epsL[0];
        for (int x = 1; x < DN_epsL.size(); x++)
        {
            if (DN_epsL[x] > max_DN)
            {
                max_DN = DN_epsL[x];
                max_index = x;
            }
        }
        // Формирование массива измеренных углов места
        eps_real[i] = eps_step * max_index + eps_start;
    }

    return eps_real;
}

// Проводка по высоте
// Дальности
QVector<double> RadarCharacteristics::Distances(double H0, const QVector<double> &eps)
{
    // H0 - Высота полета, км
    // eps - Массив углов места
    double radius = 6.371e3; // Радиус Земли (в километрах)
    double radius_eff = radius * 1.33; // Эффективный радиус Земли

    QVector<double> distance(eps.size());
    for (int i = 0; i < distance.size(); i++)
    {
        distance[i] = radius_eff *
                (-sin(eps[i] * RAD) + sqrt((pow(sin(eps[i] * RAD), 2)) + 2 * H0 / radius_eff));
    }

    return distance;
}

// Высоты
QVector<double> RadarCharacteristics::Heights(const QVector<double> &eps, const QVector<double> &eps_real,
                                              const QVector<double> & distance)
{
    // Ny - Количество элементов по x
    // dy - Расстояние между элементами по x
    // f - Частота
    // eps_a - Наклон АР
    // Hgc - Высота подъема геометрического центра АР
    // H0 - Высота полета, км
    // reflection - Отражения от земли (есть/нет)
    // eps - Массив углов места
    double radius = 6.371e3; // Радиус Земли (в километрах)
    double radius_eff = radius * 1.33; // Эффективный радиус Земли

    QVector<double> height(eps.size());
    for (int i = 0; i < height.size(); i++)
    {
        height[i] = distance[i] *
                sin(eps_real[i] * RAD) + pow(distance[i], 2) / (2 * radius_eff);
    }

    return height;
}
