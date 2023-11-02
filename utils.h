//NANOROD: utils.h Utilities Heade (Revision Date: Oct 27, 2023)
#ifndef _UTILS_H
#define _UTILS_H
#include "header.h"
#include "xyz.h"
#include "quarternion.h"
#include "system.h"
XYZ image(XYZ , double);
double min_d2(XYZ , XYZ , double );
double myfmod(double x, double y);
XYZ RandomTranslate(XYZ old, double step,double u,double v);

quarternion RandomRotate(quarternion old, double step,double a,double b,double c);
int getNum(vector<int>& v)
void generateRandom(int n)

#endif
