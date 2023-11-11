//NANOROD: utils.h Utilities Heade (Revision Date: Oct 27, 2023)
#ifndef _UTILS_H
#define _UTILS_H
#include "header.h"
#include "xyz.h"
#include "quarternion.h"
#include "system.h"
XYZ image(XYZ , double);
double min_d2(XYZ a, XYZ b, double L);
double myfmod(double x, double y);
XYZ RandomTranslate(XYZ old, double step,double u,double v);
double inner_product(XYZ a,XYZ b);
XYZ cross_product(XYZ a,XYZ b);
double angle_vectors(XYZ a,XYZ b);
double dihedral_vectors(XYZ a,XYZ b,XYZ c);
quarternion RandomRotate(quarternion old, double step,double a,double b,double c);
int getNum(vector<int>& v);
vector<int> generateRandom(int n);
int GridIndex_index(int i,int j,int k,int n);//return gridindex based on index of x, y, z (n is grids number in one direction)directions 
int GridIndex_xyz(XYZ& p,int n,double dl,double L);//return gridindex based on position p, dl is grid length
void GridLoc(int& i,int& j,int& k,int n,int index);//update xyz index based on grid index
int neighborarm(int n);
quarternion Rotate(quarternion old,double a,double b,double c);
XYZ real_vector(XYZ origin,double L);
#endif
