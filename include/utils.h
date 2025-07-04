//NANOROD: utils.h Utilities Heade (Revision Date: Oct 27, 2023)
#ifndef _UTILS_H
#define _UTILS_H
#include "header.h"
#include "xyz.h"
#include "system.h"
#include "quaternion.h"
XYZ image(XYZ , double);
double min_d2(XYZ a, XYZ b, double BoxLength);
double myfmod(double x, double y);
XYZ RandomTranslate(XYZ old, double step,double u,double v);
quaternion RandomRotation(quaternion old,double step,double u,double v);
double inner_product(XYZ a,XYZ b);
XYZ cross_product(XYZ a,XYZ b);
double angle_vectors(XYZ a,XYZ b);
double dihedral_vectors(XYZ a,XYZ b,XYZ c);
int getNum(vector<int>& v);
vector<int> generateRandom(int n);
int GridIndex_index(int i,int j,int k,int n);//return gridindex based on index of x, y, z (n is grids number in one direction)directions 
int GridIndex_xyz(XYZ& p,int n,double dl,double BoxLength);//return gridindex based on position p, dl is grid length
void GridLoc(int& i,int& j,int& k,int n,int index);//update xyz index based on grid index
int neighborarm(int n);
//quaternion Rotate(quaternion old,double a,double b,double c);
XYZ real_vector(XYZ origin,double BoxLength);
// vector<Molecule> generate_unitcell(Molecule M1);
// vector<Molecule> generate_newunit(vector<Molecule> origin_cell,int latticeindex1,int latticeindex2,int latticeindex3);

#endif
