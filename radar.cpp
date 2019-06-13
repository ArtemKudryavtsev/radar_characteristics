/****************************************************************************
**           Author: Artem Kudryavtsev                                    **
**             Date: 24.04.19                                             **
**          Version: 1.1.0                                                **
****************************************************************************/

#include "radar.h"
#include "radarcharacteristics.h"
#include "ui_radar.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <exception>

// Конструктор главного класса приложения
Radar::Radar(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::Radar)
{
    ui->setupUi(this);

    // Меню
    // Комбинация клавиш для выхода из приложения
    ui->actionExit->setShortcut(tr("Alt+Q"));
    // Установка соединения для клавиши выхода
    connect(ui->actionExit, SIGNAL(triggered()), this, SLOT(exit()));
    ui->menu->addAction(ui->actionExit);

}

// Деструктор главного класса приложения
Radar::~Radar()
{
    delete ui;
}

// Функция построения графика по двум векторам
void Radar::plot(const QVector<double> &x, const QVector<double> &y, QCustomPlot * area)
{
    // Увеличение графика
    area->setInteractions(QCP::iRangeZoom | QCP::iRangeDrag | QCP::iSelectPlottables);
    area->clearGraphs(); // Если нужно, очищаем все графики
    // Добавляем один график в widget
    area->addGraph();
    // Говорим, что отрисовать нужно график по нашим двум массивам x и y
    area->graph(0)->setData(x, y);
    // Диапазон по оси X
    area->xAxis->setRange(x[0], x[x.size()-1]);//Для оси Ox

    // Вычисление диапазона по оси Y
    double minY = y[0];
    double maxY = y[0];
    for (int i = 1; i < y.size(); i++)
    {
        if (y[i] < minY)
            minY = y[i];

        if (y[i] > maxY)
            maxY = y[i];
    }
    area->yAxis->setRange(minY, maxY);
    // И перерисуем график на нашем widget
    area->replot();
}

/*
 * ДН ПО АЗИМУТУ
 */

// Параметры для ДН по азимуту по умолчанию
void Radar::on_buttonDefaultBeta_clicked()
{
    ui->inputNxBeta->setText("20"); // Количество столбцов
    ui->inputDxBeta->setText("1.23"); // Расстояние между элементами
    ui->inputFBeta->setText("123"); // Частота
    ui->inputTargetBeta->setText("0.0"); // Угол цели
}

// Построение ДН по азимуту
void Radar::on_buttonPlotBeta_clicked()
{
    Nx = ui->inputNxBeta->text().toDouble(); // Количество столбцов
    dx = ui->inputDxBeta->text().toDouble(); // Расстояние между элементами
    f_beta = ui->inputFBeta->text().toDouble() * 1e6; // Частота
    beta_target = ui->inputTargetBeta->text().toDouble(); // Угол цели

    double beta_start = -30.0; // Начало сканирования
    double beta_step = 0.1; // Шаг сканирования
    double beta_stop = 30.0; // Конец сканирования

    int beta_size = (beta_stop - beta_start) / beta_step + 1; // размер массива углов

    beta.resize(beta_size); // Вектор азимутов
    beta[0] = beta_start;
    for (int i = 1; i < beta.size(); i++)
    {
        beta[i] = beta[i-1] + beta_step;
    }

    DN_beta.resize(beta_size); // ДН (линейный масштаб)

    RadarCharacteristics radar; // Создаем объект класса RadarCharacteristics,
                                // чтобы получить доступ к методу DN_beta
    DN_beta = radar.DN_beta(Nx, dx, f_beta, beta_target, beta);

    double DN_beta_max = DN_beta[0]; // Максимальное значение ДН
    for (QVector<double>::iterator p = DN_beta.begin(); p != DN_beta.end(); p++)
    {
        if (*p > DN_beta_max)
            DN_beta_max = *p;
    }

    // Нормируем массив ДН к единице
    for (auto p = DN_beta.begin(); p != DN_beta.end(); p++)
    {
        *p = *p / DN_beta_max;
    }

    // Если логарифмический масштаб, то преобразуем ДН к логарифмическому масштабу
    if (ui->radioButtonLogBeta->isChecked())
    {
        for (auto p = DN_beta.begin(); p != DN_beta.end(); p++)
        {
            *p = 20 * log10(*p);
        }
    }
    QCustomPlot * area = ui->plotAreaBeta; // Объект графика
    // Построение графика
    // Подписываем оси Ox и Oy
    if (ui->radioButtonLinBeta->isChecked())
        area->yAxis->setLabel("Амплитуда");
    else
        area->yAxis->setLabel("Амплитуда, дБ");
    area->xAxis->setLabel("Азимут, градусы");
    this->plot(beta, DN_beta, area); 
}

// Сохранение в файл
void Radar::on_buttonSaveBeta_clicked()
{
    // Запись в файл через диалоговое окно с выбором пути
    QString name;
    if (ui->radioButtonLinBeta->isChecked())
        name = "DN_beta_lin.txt";
    else
        name = "DN_beta_log.txt";

    std::string str = QFileDialog::getSaveFileName(0, tr("Сохранить в файл"), name, tr("*.txt")).toStdString();
    // Запись в стиле C++
    std::ofstream fout;
    fout.open(str);
    if (fout.is_open()) {

        // Установка фиксированной формы вывода
        fout.setf(std::ios::fixed);
        if (ui->radioButtonLinBeta->isChecked())
            fout << "Азимут" << '\t' << "Амплитуда\n";
        else
            fout << "Азимут" << '\t' << "Амплитуда, дБ\n";
        for (int i = 0; i < beta.size(); i++)
            fout << std::setprecision(1) << beta[i] << "\t" << std::setprecision(6) << DN_beta[i] << '\n';

    }
    fout.close();
}

/*
 * ДН ПО УГЛУ МЕСТА
 */

// Параметры для ДН по углу места по умолчанию
void Radar::on_buttonDefaultEps_clicked()
{
    ui->inputNyEps->setText("10"); // Количество строк
    ui->inputDyEps->setText("1.23"); // Расстояние между элементами
    ui->inputFEps->setText("123"); // Частота
    ui->inputAngleEps->setText("0"); // Наклон АР
    ui->inputHgcEps->setText("10"); // Высота поднятия геометрического центра АР
    ui->inputTargetEps->setText("0.1"); // Угол цели
}

// Построение ДН по углу места
void Radar::on_buttonPlotEps_clicked()
{
    Ny = ui->inputNyEps->text().toDouble(); // Количество строк
    dy = ui->inputDyEps->text().toDouble(); // Расстояние между элементами
    f_eps = ui->inputFEps->text().toDouble() * 1e6; // Частота
    eps_a = ui->inputAngleEps->text().toDouble(); // Наклон АР
    Hgc = ui->inputHgcEps->text().toDouble(); // Высота подъема геометрического центра
    eps_target = ui->inputTargetEps->text().toDouble(); // Угол цели

    bool reflection = false; // Нет отражений от земли
    if (ui->checkRefractionEps->isChecked())
        reflection = true; // Отражения есть

    double eps_start = 0.1; // Начало сканирования
    double eps_step = 0.1; // Шаг сканирования
    double eps_stop = 30.0; // Конец сканирования

    int eps_size = (eps_stop - eps_start) / eps_step + 1; // размер массива углов

    eps.resize(eps_size); // Вектор углов места
    eps[0] = eps_start;
    for (int i = 1; i < eps.size(); i++)
    {
        eps[i] = eps[i-1] + eps_step;
    }

    RadarCharacteristics radar; // Создаем объект класса RadarCharacteristics,
    // чтобы получить доступ к методу DN_beta
    DN_eps.resize(eps_size); // ДН (линейный масштаб)
    DN_eps = radar.DN_eps(Ny, dy, f_eps, eps_a, eps_target, Hgc, reflection, eps);

    double DN_eps_max = DN_eps[0]; // Максимальное значение ДН
    for (QVector<double>::iterator p = DN_eps.begin(); p != DN_eps.end(); p++)
    {
        if (*p > DN_eps_max)
            DN_eps_max = *p;
    }

    // Нормируем массив ДН к единице
    for (auto p = DN_eps.begin(); p != DN_eps.end(); p++)
    {
        *p = *p / DN_eps_max;
    }

    // Если логарифмический масштаб, то преобразуем ДН к логарифмическому масштабу
    if (ui->radioButtonLogEps->isChecked())
    {
        for (auto p = DN_eps.begin(); p != DN_eps.end(); p++)
        {
            *p = 20 * log10(*p);
        }
    }

    QCustomPlot * area = ui->plotAreaEps; // Объект графика
    // Построение графика
    // Подписываем оси Ox и Oy
    if (ui->radioButtonLinEps->isChecked())
        area->yAxis->setLabel("Амплитуда");
    else
        area->yAxis->setLabel("Амплитуда, дБ");
    area->xAxis->setLabel("Угол места, градусы");
    this->plot(eps, DN_eps, area);
}

// Сохранение в файл
void Radar::on_buttonSaveEps_clicked()
{
    // Запись в файл через диалоговое окно с выбором пути
    QString name;
    if (ui->radioButtonLinEps->isChecked() && ui->checkRefractionEps->isChecked())
        name = "DN_eps_with_earth_lin.txt";
    else if (ui->radioButtonLinEps->isChecked() && !ui->checkRefractionEps->isChecked())
        name = "DN_eps_without_earth_lin.txt";
    else if (ui->radioButtonLogEps->isChecked() && ui->checkRefractionEps->isChecked())
        name = "DN_eps_with_earth_log.txt";
    else if (ui->radioButtonLogEps->isChecked() && !ui->checkRefractionEps->isChecked())
        name = "DN_eps_without_earth_log.txt";

    // Запись в файл через диалоговое окно с выбором пути
    std::string str = QFileDialog::getSaveFileName(0, tr("Сохранить в файл"), name, tr("*.txt")).toStdString();
    // Запись в стиле C++
    std::ofstream fout;
    fout.open(str);
    if (fout.is_open()) {

        // Установка фиксированной формы вывода
        fout.setf(std::ios::fixed);
        if (ui->radioButtonLinEps->isChecked())
            fout << "УМ" << '\t' << "Амплитуда\n";
        else
            fout << "УМ" << '\t' << "Амплитуда, дБ\n";
        for (int i = 0; i < eps.size(); i++)
            fout << std::setprecision(1) << eps[i] << '\t' << std::setprecision(6) << DN_eps[i] << '\n';

    }
    fout.close();
}

/*
 * ПЕЛЕНГАЦИОННАЯ ХАРАКТЕРИСТИКА
 */

// Параметры для пеленгационной характеристики по умолчанию
void Radar::on_buttonDefaultPeleng_clicked()
{
    ui->inputNyPeleng->setText("10"); // Количество строк
    ui->inputDyPeleng->setText("1.23"); // Расстояние между элементами
    ui->inputFPeleng->setText("123"); // Частота
    ui->inputAnglePeleng->setText("20"); // Наклон АР
    ui->inputHgcPeleng->setText("10"); // Высота поднятия геометрического центра АР
}

// Построение пеленгационной характеристики
void Radar::on_buttonPlotPeleng_clicked()
{
    Ny_pel = ui->inputNyPeleng->text().toDouble(); // Количество строк
    dy_pel = ui->inputDyPeleng->text().toDouble(); // Расстояние между элементами
    f_eps_pel = ui->inputFPeleng->text().toDouble() * 1e6; // Частота
    eps_a_pel = ui->inputAnglePeleng->text().toDouble(); // Наклон АР
    Hgc_pel = ui->inputHgcPeleng->text().toDouble(); // Высота подъема геометрического центра

    bool reflection = false; // Отражений от земли нет
    if (ui->checkRefractionPeleng->isChecked())
        reflection = true; // Отражения от земли есть

    RadarCharacteristics radar; // Создаем объект класса RadarCharacteristics,
    // чтобы получить доступ к методу DN_beta

    double eps_start = 0.1; // Начальный угол сканирования
    double eps_step = 0.1; // Шаг сканирования
    double eps_stop = 70.0; // Конечный угол сканирования
    int eps_size = (eps_stop - eps_start) / eps_step + 1;

    eps_orig_pel.resize(eps_size); // Вектор углов места
    eps_orig_pel[0] = eps_start;
    for (int i = 1; i < eps_orig_pel.size(); i++)
    {
        eps_orig_pel[i] = eps_orig_pel[i-1] + eps_step;
    }

    eps_real_pel.resize(eps_size); // Измеренные углы места
    eps_real_pel = radar.Peleng(Ny_pel, dy_pel, f_eps_pel, eps_a_pel, Hgc_pel, reflection, eps_orig_pel);

    // Увеличение графика
    ui->plotAreaPeleng->setInteractions(QCP::iRangeZoom | QCP::iRangeDrag | QCP::iSelectPlottables);
    ui->plotAreaPeleng->clearGraphs(); // Если нужно, очищаем все графики
    // Добавляем один график в widget
    ui->plotAreaPeleng->addGraph();
    // Говорим, что отрисовать нужно график по нашим двум массивам x и y
    ui->plotAreaPeleng->graph(0)->setData(eps_orig_pel, eps_real_pel);

    ui->plotAreaPeleng->xAxis->setLabel("Угол места истинный, градус");
    ui->plotAreaPeleng->yAxis->setLabel("Угол места измеренный, градус");

    // Диапазон по оси X
    ui->plotAreaPeleng->xAxis->setRange(eps_orig_pel[0], eps_orig_pel[eps_orig_pel.size() - 1]);//Для оси Ox
    ui->plotAreaPeleng->yAxis->setRange(eps_real_pel[0], eps_real_pel[eps_real_pel.size() - 1]);
    // И перерисуем график на нашем widget
    ui->plotAreaPeleng->replot();
}

// Сохранение в файл
void Radar::on_buttonSavePeleng_clicked()
{
    // Запись в файл через диалоговое окно с выбором пути
    QString name;
    if (ui->checkRefractionPeleng->isChecked())
        name = "Peleng_with_earth.txt";
    else
        name = "Peleng_without_earth.txt";

    // Запись в файл через диалоговое окно с выбором пути
    std::string str = QFileDialog::getSaveFileName(0, tr("Сохранить в файл"), name, tr("*.txt")).toStdString();
    // Запись в стиле C++
    std::ofstream fout;
    fout.open(str);
    if (fout.is_open()) {

        // Установка фиксированной формы вывода
        fout.setf(std::ios::fixed);
        fout << "УМ истинный" << '\t' << "УМ измеренный\n";
        for (int i = 0; i < eps_orig_pel.size(); i++)
            fout << std::setprecision(2) << eps_orig_pel[i] << "\t\t" << std::setprecision(2) << eps_real_pel[i] << '\n';

    }
    fout.close();
}

/*
 * ПРОВОДКА ПО ВЫСОТЕ
 */

// Параметры для проводки по высоте по умолчанию
void Radar::on_buttonDefaultProvodka_clicked()
{
    ui->inputNyProvodka->setText("10"); // Количество строк
    ui->inputDyProvodka->setText("1.23"); // Расстояние между элементами
    ui->inputFProvodka->setText("123"); // Частота
    ui->inputAngleProvodka->setText("20"); // Наклон АР
    ui->inputHgcProvodka->setText("10"); // Высота поднятия геометрического центра АР
    ui->inputH0Provodka->setText("10"); // Высота полета
}

// Построение проводки по высоте
void Radar::on_buttonPlotProvodka_clicked()
{
    Ny_prov = ui->inputNyProvodka->text().toDouble(); // Количество строк
    dy_prov = ui->inputDyProvodka->text().toDouble(); // Расстояние между элементами
    f_eps_prov = ui->inputFProvodka->text().toDouble() * 1e6; // Частота
    eps_a_prov = ui->inputAngleProvodka->text().toDouble(); // Наклон АР
    Hgc_prov = ui->inputHgcProvodka->text().toDouble(); // Высота подъема геометрического центра
    H0_prov = ui->inputH0Provodka->text().toDouble(); // Высота полета

    bool reflection = false; // Нет отражений от земли
    if (ui->checkRefractionProvodka->isChecked())
        reflection = true; // Отражения есть

    RadarCharacteristics radar; // Создаем объект класса RadarCharacteristics,
    // чтобы получить доступ к методу DN_beta

    double eps_start = 0.1; // Начальный угол сканирования
    double eps_step = 0.1; // Шаг сканирования
    double eps_stop = 70.0; // Конечный угол сканирования
    int eps_size = (eps_stop - eps_start) / eps_step + 1;

    eps_orig_prov.resize(eps_size); // Вектор углов места
    eps_orig_prov[0] = eps_start;
    for (int i = 1; i < eps_orig_prov.size(); i++)
    {
        eps_orig_prov[i] = eps_orig_prov[i-1] + eps_step;
    }

    eps_real_prov = radar.Peleng(Ny_prov, dy_prov, f_eps_prov, eps_a_prov, Hgc_prov, reflection, eps_orig_prov);

    // Дальности
    distance_prov = radar.Distances(H0_prov, eps_orig_prov);
    // Высоты
    height_prov = radar.Heights(eps_orig_prov, eps_real_prov, distance_prov);

    // Увеличение графика
    ui->plotAreaProvodka->setInteractions(QCP::iRangeZoom | QCP::iRangeDrag | QCP::iSelectPlottables);
    ui->plotAreaProvodka->clearGraphs();//Если нужно, но очищаем все графики
    //Добавляем один график в widget
    ui->plotAreaProvodka->addGraph();
    //Говорим, что отрисовать нужно график по нашим двум массивам x и y
    ui->plotAreaProvodka->graph(0)->setData(distance_prov, height_prov);

    //Подписываем оси Ox и Oy
    ui->plotAreaProvodka->xAxis->setLabel("Дальность, км");
    ui->plotAreaProvodka->yAxis->setLabel("Высота, км");

    //Установим область, которая будет показываться на графике
    ui->plotAreaProvodka->xAxis->setRange(distance_prov[0], distance_prov[distance_prov.size() - 1]);//Для оси Ox

    double minY = height_prov[0], maxY = height_prov[0];
    for (int i = 1; i < height_prov.size(); i++)
    {
        if (height_prov[i]<minY) minY = height_prov[i];
        if (height_prov[i]>maxY) maxY = height_prov[i];
    }
    ui->plotAreaProvodka->yAxis->setRange(minY, maxY);//Для оси Oy
    //И перерисуем график на нашем widget
    ui->plotAreaProvodka->replot();
}

// Сохранение в файл
void Radar::on_buttonSaveProvodka_clicked()
{
    // Запись в файл через диалоговое окно с выбором пути
    QString name;
    if (ui->checkRefractionProvodka->isChecked())
        name = "Provodka_with_earth.txt";
    else
        name = "Provodka_without_earth.txt";

    // Запись в файл через диалоговое окно с выбором пути
    std::string str = QFileDialog::getSaveFileName(0, tr("Сохранить в файл"), name, tr("*.txt")).toStdString();
    // Запись в стиле C++
    std::ofstream fout;
    fout.open(str);
    if (fout.is_open()) {

        // Установка фиксированной формы вывода
        fout.setf(std::ios::fixed);
        fout << "Дальность, км" << '\t' << "Высота, км\n";
        for (int i = 0; i < distance_prov.size(); i++)
            fout << std::setprecision(2) << distance_prov[i] << "\t\t" << std::setprecision(2) << height_prov[i] << '\n';

    }
    fout.close();
}
