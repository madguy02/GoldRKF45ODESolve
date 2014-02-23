/* Copyright (c) <2014> Author Vance King Saxbe. A, and contributors Power Dominion Enterprise, Precieux Consulting and other contributors. Modelled, Architected and designed by Vance King Saxbe. A. with the geeks from GoldSax Consulting and GoldSax Technologies email @vsaxbe@yahoo.com. Development teams from Power Dominion Enterprise, Precieux Consulting. Project sponsored by GoldSax Foundation, GoldSax Group and executed by GoldSax Manager.*//* Copyright (c) <2014> Author Vance King Saxbe. A, and contributors Power Dominion Enterprise, Precieux Consulting and other contributors. Modelled, Architected and designed by Vance King Saxbe. A. with the geeks from GoldSax Consulting and GoldSax Technologies email @vsaxbe@yahoo.com. Development teams from Power Dominion Enterprise, Precieux Consulting. Project sponsored by GoldSax Foundation, GoldSax Group and executed by GoldSax Manager.*/float r4_abs ( float x );
float r4_epsilon ( );
void r4_fehl ( void f ( float t, float y[], float yp[] ), int neqn, 
  float y[], float t, float h, float yp[], float f1[], float f2[], float f3[], 
  float f4[], float f5[], float s[] );
float r4_max ( float x, float y );
float r4_min ( float x, float y );
int r4_GoldSaxRKF45ODESolve ( void f ( float t, float y[], float yp[] ), int neqn,
  float y[], float yp[], float *t, float tout, float *relerr, float abserr, 
  int flag );
float r4_sign ( float x );

double r8_abs ( double x );
double r8_epsilon ( );
void r8_fehl ( void f ( double t, double y[], double yp[] ), int neqn, 
  double y[], double t, double h, double yp[], double f1[], double f2[], double f3[], 
  double f4[], double f5[], double s[] );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
int r8_GoldSaxRKF45ODESolve ( void f ( double t, double y[], double yp[] ), int neqn,
  double y[], double yp[], double *t, double tout, double *relerr, double abserr, 
  int flag );
double r8_sign ( double x );

void timestamp ( );
