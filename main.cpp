#include "radar.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Radar w;
    w.show();

    return a.exec();
}
