/* Copyright (c) <2014> Author Vance King Saxbe. A, and contributors Power Dominion Enterprise, Precieux Consulting and other contributors. Modelled, Architected and designed by Vance King Saxbe. A. with the geeks from GoldSax Consulting and GoldSax Technologies email @vsaxbe@yahoo.com. Development teams from Power Dominion Enterprise, Precieux Consulting. Project sponsored by GoldSax Foundation, GoldSax Group and executed by GoldSax Manager.*/# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "GoldSaxRKF45ODESolve.hpp"


float r4_abs ( float x )


{
  if ( 0.0 <= x )
  {
    return x;
  } 
  else
  {
    return ( -x );
  }
}

float r4_epsilon ( )


{
  float r;

  r = 1.0;

  while ( 1.0 < ( float ) ( 1.0 + r )  )
  {
    r = r / 2.0;
  }

  return ( 2.0 * r );
}

void r4_fehl ( void f ( float t, float y[], float yp[] ), int neqn, 
  float y[], float t, float h, float yp[], float f1[], float f2[], float f3[], 
  float f4[], float f5[], float s[] )


{
  float ch;
  int i;

  ch = h / 4.0;

  for ( i = 0; i < neqn; i++ )
  {
    f5[i] = y[i] + ch * yp[i];
  }

  f ( t + ch, f5, f1 );

  ch = 3.0 * h / 32.0;

  for ( i = 0; i < neqn; i++ )
  {
    f5[i] = y[i] + ch * ( yp[i] + 3.0 * f1[i] );
  }

  f ( t + 3.0 * h / 8.0, f5, f2 );

  ch = h / 2197.0;

  for ( i = 0; i < neqn; i++ )
  {
    f5[i] = y[i] + ch * 
    ( 1932.0 * yp[i] 
    + ( 7296.0 * f2[i] - 7200.0 * f1[i] ) 
    );
  }

  f ( t + 12.0 * h / 13.0, f5, f3 );

  ch = h / 4104.0;

  for ( i = 0; i < neqn; i++ )
  {
    f5[i] = y[i] + ch * 
    ( 
      ( 8341.0 * yp[i] - 845.0 * f3[i] ) 
    + ( 29440.0 * f2[i] - 32832.0 * f1[i] ) 
    );
  }

  f ( t + h, f5, f4 );

  ch = h / 20520.0;

  for ( i = 0; i < neqn; i++ )
  {
    f1[i] = y[i] + ch * 
    ( 
      ( -6080.0 * yp[i] 
      + ( 9295.0 * f3[i] - 5643.0 * f4[i] ) 
      ) 
    + ( 41040.0 * f1[i] - 28352.0 * f2[i] ) 
    );
  }

  f ( t + h / 2.0, f1, f5 );

  ch = h / 7618050.0;

  for ( i = 0; i < neqn; i++ )
  {
    s[i] = y[i] + ch * 
    ( 
      ( 902880.0 * yp[i] 
      + ( 3855735.0 * f3[i] - 1371249.0 * f4[i] ) ) 
    + ( 3953664.0 * f2[i] + 277020.0 * f5[i] ) 
    );
  }

  return;
}

float r4_max ( float x, float y )


{
  if ( y < x )
  {
    return x;
  } 
  else
  {
    return y;
  }
}

float r4_min ( float x, float y )

{
  if ( y < x )
  {
    return y;
  } 
  else
  {
    return x;
  }
}

int r4_GoldSaxRKF45ODESolve ( void f ( float t, float y[], float yp[] ), int neqn,
  float y[], float yp[], float *t, float tout, float *relerr, float abserr, 
  int flag )

{
# define MAXNFE 3000

  static float abserr_save = -1.0;
  float ae;
  float dt;
  float ee;
  float eeoet;
  float eps;
  float esttol;
  float et;
  float *f1;
  float *f2;
  float *f3;
  float *f4;
  float *f5;
  int flag_return;
  static int flag_save = -1000;
  static float h = -1.0;
  bool hfaild;
  float hmin;
  int i;
  static int init = -1000;
  int k;
  static int kflag = -1000;
  static int kop = -1;
  int mflag;
  static int nfe = -1;
  bool output;
  float relerr_min;
  static float relerr_save = -1.0;
  static float remin = 1.0E-12;
  float s;
  float scale;
  float tol;
  float toln;
  float ypk;

  flag_return = flag;

  eps = r4_epsilon ( );

  if ( neqn < 1 )
  {
    flag_return = 8;
    return flag_return;
  }

  if ( (*relerr) < 0.0 )
  {
    flag_return = 8;
    return flag_return;
  }

  if ( abserr < 0.0 )
  {
    flag_return = 8;
    return flag_return;
  }

  if ( flag_return == 0 || 8 < flag_return  || flag_return < -2 )
  {
    flag_return = 8;
    return flag_return;
  }

  mflag = abs ( flag_return );

  if ( mflag != 1 )
  {
    if ( *t == tout && kflag != 3 )
    {
      flag_return = 8;
      return flag_return;
    }
//
//  FLAG = -2 or +2:
//
    if ( mflag == 2 )
    {
      if ( kflag == 3 )
      {
        flag_return = flag_save;
        mflag = abs ( flag_return );
      }
      else if ( init == 0 )
      {
        flag_return = flag_save;
      }
      else if ( kflag == 4 )
      {
        nfe = 0;
      }
      else if ( kflag == 5 && abserr == 0.0 )
      {
        cerr << "\n";
        cerr << "R4_GoldSaxRKF45ODESolve - Fatal error!\n";
        cerr << "  KFLAG = 5 and ABSERR = 0.0\n";
        exit ( 1 );
      }
      else if ( kflag == 6 && (*relerr) <= relerr_save && abserr <= abserr_save )
      {
        cerr << "\n";
        cerr << "R4_GoldSaxRKF45ODESolve - Fatal error!\n";
        cerr << "  KFLAG = 6 and\n";
        cerr << "  RELERR <= RELERR_SAVE and\n";
        cerr << "  ABSERR <= ABSERR_SAVE\n";
        exit ( 1 );
      }
    }
//
//  FLAG = 3, 4, 5, 6, 7 or 8.
//
    else
    {
      if ( flag_return == 3 )
      {
        flag_return = flag_save;
        if ( kflag == 3 )
        {
          mflag = abs ( flag_return );
        }
      }
      else if ( flag_return == 4 )
      {
        nfe = 0;
        flag_return = flag_save;
        if ( kflag == 3 )
        {
          mflag = abs ( flag_return );
        }
      }
      else if ( flag_return == 5 && 0.0 < abserr )
      {
        flag_return = flag_save;
        if ( kflag == 3 )
        {
          mflag = abs ( flag_return );
        }
      }
//
//  Integration cannot be continued because the user did not respond to
//  the instructions pertaining to FLAG = 5, 6, 7 or 8.
//
      else
      {
        cerr << "\n";
        cerr << "R4_GoldSaxRKF45ODESolve - Fatal error!\n";
        cerr << "  Integration cannot be continued.\n";
        cerr << "  The user did not respond to the output\n";
        cerr << "  value FLAG = 5, 6, 7, or 8.\n";
        exit ( 1 );
      }
    }
  }

  flag_save = flag_return;
  kflag = 0;
//
//  Save RELERR and ABSERR for checking input on subsequent calls.
//
  relerr_save = (*relerr);
  abserr_save = abserr;
//
//  Restrict the relative error tolerance to be at least 
//
//    2*EPS+REMIN 
//
//  to avoid limiting precision difficulties arising from impossible 
//  accuracy requests.
//
  relerr_min = 2.0 * r4_epsilon ( ) + remin;
//
//  Is the relative error tolerance too small?
//
  if ( (*relerr) < relerr_min )
  {
    (*relerr) = relerr_min;
    kflag = 3;
    flag_return = 3;
    return flag_return;
  }

  dt = tout - *t;

  f1 = new float[neqn];
  f2 = new float[neqn];
  f3 = new float[neqn];
  f4 = new float[neqn];
  f5 = new float[neqn];

  if ( mflag == 1 )
  {
    init = 0;
    kop = 0;
    f ( *t, y, yp );
    nfe = 1;

    if ( *t == tout )
    {
      flag_return = 2;
      return flag_return;
    }

  }

  if ( init == 0 )
  {
    init = 1;
    h = r4_abs ( dt );
    toln = 0.0;

    for ( k = 0; k < neqn; k++ )
    {
      tol = (*relerr) * r4_abs ( y[k] ) + abserr;
      if ( 0.0 < tol )
      {
        toln = tol;
        ypk = r4_abs ( yp[k] );
        if ( tol < ypk * pow ( h, 5 ) )
        {
          h = ( float ) pow ( ( double ) ( tol / ypk ), 0.2 );
        }
      }
    }

    if ( toln <= 0.0 )
    {
      h = 0.0;
    }

    h = r4_max ( h, 26.0 * eps * r4_max ( r4_abs ( *t ), r4_abs ( dt ) ) );

    if ( flag_return < 0 )
    {
      flag_save = -2;
    }
    else
    {
      flag_save = 2;
    }
  }
//
//  Set stepsize for integration in the direction from T to TOUT.
//
  h = r4_sign ( dt ) * r4_abs ( h );
//
//  Test to see if too may output points are being requested.
//
  if ( 2.0 * r4_abs ( dt ) <= r4_abs ( h ) )
  {
    kop = kop + 1;
  }
//
//  Unnecessary frequency of output.
//
  if ( kop == 100 )
  {
    kop = 0;
    delete [] f1;
    delete [] f2;
    delete [] f3;
    delete [] f4;
    delete [] f5;
    flag_return = 7;
    return flag_return;
  }
//
//  If we are too close to the output point, then simply extrapolate and return.
//
  if ( r4_abs ( dt ) <= 26.0 * eps * r4_abs ( *t ) )
  {
    *t = tout;
    for ( i = 0; i < neqn; i++ )
    {
      y[i] = y[i] + dt * yp[i];
    }
    f ( *t, y, yp );
    nfe = nfe + 1;

    delete [] f1;
    delete [] f2;
    delete [] f3;
    delete [] f4;
    delete [] f5;
    flag_return = 2;
    return flag_return;
  }
//
//  Initialize the output point indicator.
//
  output = false;
//
//  To avoid premature underflow in the error tolerance function,
//  scale the error tolerances.
//
  scale = 2.0 / (*relerr);
  ae = scale * abserr;
//
//  Step by step integration.
//
  for ( ; ; )
  {
    hfaild = false;
//
//  Set the smallest allowable stepsize.
//
    hmin = 26.0 * eps * r4_abs ( *t );

    dt = tout - *t;

    if ( 2.0 * r4_abs ( h ) <= r4_abs ( dt ) )
    {
    }
    else
//
//  Will the next successful step complete the integration to the output point?
//
    {
      if ( r4_abs ( dt ) <= r4_abs ( h ) )
      {
        output = true;
        h = dt;
      }
      else
      {
        h = 0.5 * dt;
      }

    }

    for ( ; ; )
    {

      if ( MAXNFE < nfe )
      {
        kflag = 4;
        delete [] f1;
        delete [] f2;
        delete [] f3;
        delete [] f4;
        delete [] f5;
        flag_return = 4;
        return flag_return;
      }

      r4_fehl ( f, neqn, y, *t, h, yp, f1, f2, f3, f4, f5, f1 );
      nfe = nfe + 5;

      eeoet = 0.0;
 
      for ( k = 0; k < neqn; k++ )
      {
        et = r4_abs ( y[k] ) + r4_abs ( f1[k] ) + ae;

        if ( et <= 0.0 )
        {
          delete [] f1;
          delete [] f2;
          delete [] f3;
          delete [] f4;
          delete [] f5;
          flag_return = 5;
          return flag_return;
        }

        ee = r4_abs 
        ( ( -2090.0 * yp[k] 
          + ( 21970.0 * f3[k] - 15048.0 * f4[k] ) 
          ) 
        + ( 22528.0 * f2[k] - 27360.0 * f5[k] ) 
        );

        eeoet = r4_max ( eeoet, ee / et );

      }

      esttol = r4_abs ( h ) * eeoet * scale / 752400.0;

      if ( esttol <= 1.0 )
      {
        break;
      }

      hfaild = true;
      output = false;

      if ( esttol < 59049.0 )
      {
        s = 0.9 / ( float ) pow ( ( double ) esttol, 0.2 );
      }
      else
      {
        s = 0.1;
      }

      h = s * h;

      if ( r4_abs ( h ) < hmin )
      {
        kflag = 6;
        delete [] f1;
        delete [] f2;
        delete [] f3;
        delete [] f4;
        delete [] f5;
        flag_return = 6;
        return flag_return;
      }

    }

    *t = *t + h;
    for ( i = 0; i < neqn; i++ )
    {
      y[i] = f1[i];
    }
    f ( *t, y, yp );
    nfe = nfe + 1;

    if ( 0.0001889568 < esttol )
    {
      s = 0.9 / ( float ) pow ( ( double ) esttol, 0.2 );
    }
    else
    {
      s = 5.0;
    }

    if ( hfaild )
    {
      s = r4_min ( s, 1.0 );
    }

    h = r4_sign ( h ) * r4_max ( s * r4_abs ( h ), hmin );

    if ( output )
    {
      *t = tout;
      delete [] f1;
      delete [] f2;
      delete [] f3;
      delete [] f4;
      delete [] f5;
      flag_return = 2;
      return flag_return;
    }

    if ( flag_return <= 0 )
    {
      delete [] f1;
      delete [] f2;
      delete [] f3;
      delete [] f4;
      delete [] f5;
      flag_return = -2;
      return flag_return;
    }

  }
# undef MAXNFE
}

float r4_sign ( float x )


{
  if ( x < 0.0 )
  {
    return ( -1.0 );
  } 
  else
  {
    return ( +1.0 );
  }
}

double r8_abs ( double x )

{
  if ( 0.0 <= x )
  {
    return x;
  } 
  else
  {
    return ( -x );
  }
}

double r8_epsilon ( )


{
  double r;

  r = 1.0;

  while ( 1.0 < ( double ) ( 1.0 + r )  )
  {
    r = r / 2.0;
  }

  return ( 2.0 * r );
}

void r8_fehl ( void f ( double t, double y[], double yp[] ), int neqn, 
  double y[], double t, double h, double yp[], double f1[], double f2[], 
  double f3[], double f4[], double f5[], double s[] )

{
  double ch;
  int i;

  ch = h / 4.0;

  for ( i = 0; i < neqn; i++ )
  {
    f5[i] = y[i] + ch * yp[i];
  }

  f ( t + ch, f5, f1 );

  ch = 3.0 * h / 32.0;

  for ( i = 0; i < neqn; i++ )
  {
    f5[i] = y[i] + ch * ( yp[i] + 3.0 * f1[i] );
  }

  f ( t + 3.0 * h / 8.0, f5, f2 );

  ch = h / 2197.0;

  for ( i = 0; i < neqn; i++ )
  {
    f5[i] = y[i] + ch * 
    ( 1932.0 * yp[i] 
    + ( 7296.0 * f2[i] - 7200.0 * f1[i] ) 
    );
  }

  f ( t + 12.0 * h / 13.0, f5, f3 );

  ch = h / 4104.0;

  for ( i = 0; i < neqn; i++ )
  {
    f5[i] = y[i] + ch * 
    ( 
      ( 8341.0 * yp[i] - 845.0 * f3[i] ) 
    + ( 29440.0 * f2[i] - 32832.0 * f1[i] ) 
    );
  }

  f ( t + h, f5, f4 );

  ch = h / 20520.0;

  for ( i = 0; i < neqn; i++ )
  {
    f1[i] = y[i] + ch * 
    ( 
      ( -6080.0 * yp[i] 
      + ( 9295.0 * f3[i] - 5643.0 * f4[i] ) 
      ) 
    + ( 41040.0 * f1[i] - 28352.0 * f2[i] ) 
    );
  }

  f ( t + h / 2.0, f1, f5 );

  ch = h / 7618050.0;

  for ( i = 0; i < neqn; i++ )
  {
    s[i] = y[i] + ch * 
    ( 
      ( 902880.0 * yp[i] 
      + ( 3855735.0 * f3[i] - 1371249.0 * f4[i] ) ) 
    + ( 3953664.0 * f2[i] + 277020.0 * f5[i] ) 
    );
  }

  return;
}


double r8_max ( double x, double y )

{
  if ( y < x )
  {
    return x;
  } 
  else
  {
    return y;
  }
}


double r8_min ( double x, double y )

{
  if ( y < x )
  {
    return y;
  } 
  else
  {
    return x;
  }
}

int r8_GoldSaxRKF45ODESolve ( void f ( double t, double y[], double yp[] ), int neqn,
  double y[], double yp[], double *t, double tout, double *relerr, 
  double abserr, int flag )

{
# define MAXNFE 3000

  static double abserr_save = -1.0;
  double ae;
  double dt;
  double ee;
  double eeoet;
  double eps;
  double esttol;
  double et;
  double *f1;
  double *f2;
  double *f3;
  double *f4;
  double *f5;
  int flag_return;
  static int flag_save = -1000;
  static double h = -1.0;
  bool hfaild;
  double hmin;
  int i;
  static int init = -1000;
  int k;
  static int kflag = -1000;
  static int kop = -1;
  int mflag;
  static int nfe = -1;
  bool output;
  double relerr_min;
  static double relerr_save = -1.0;
  static double remin = 1.0E-12;
  double s;
  double scale;
  double tol;
  double toln;
  double ypk;

  flag_return = flag;

  eps = r8_epsilon ( );

  if ( neqn < 1 )
  {
    flag_return = 8;
    cerr << "\n";
    cerr << "R8_GoldSaxRKF45ODESolve - Fatal error!\n";
    cerr << "  Invalid input value of NEQN.\n";
    return flag_return;
  }

  if ( (*relerr) < 0.0 )
  {
    flag_return = 8;
    cerr << "\n";
    cerr << "R8_GoldSaxRKF45ODESolve - Fatal error!\n";
    cerr << "  Invalid input value of RELERR.\n";
    return flag_return;
  }

  if ( abserr < 0.0 )
  {
    flag_return = 8;
    cerr << "\n";
    cerr << "R8_GoldSaxRKF45ODESolve - Fatal error!\n";
    cerr << "  Invalid input value of ABSERR.\n";
    return flag_return;
  }

  if ( flag_return == 0 || 8 < flag_return  || flag_return < -2 )
  {
    flag_return = 8;
    cerr << "\n";
    cerr << "R8_GoldSaxRKF45ODESolve - Fatal error!\n";
    cerr << "  Invalid input.\n";
    return flag_return;
  }

  mflag = abs ( flag_return );
//
//  Is this a continuation call?
//
  if ( mflag != 1 )
  {
    if ( *t == tout && kflag != 3 )
    {
      flag_return = 8;
      return flag_return;
    }
//
//  FLAG = -2 or +2:
//
    if ( mflag == 2 )
    {
      if ( kflag == 3 )
      {
        flag_return = flag_save;
        mflag = abs ( flag_return );
      }
      else if ( init == 0 )
      {
        flag_return = flag_save;
      }
      else if ( kflag == 4 )
      {
        nfe = 0;
      }
      else if ( kflag == 5 && abserr == 0.0 )
      {
        cerr << "\n";
        cerr << "R8_GoldSaxRKF45ODESolve - Fatal error!\n";
        cerr << "  KFLAG = 5 and ABSERR = 0.0\n";
        exit ( 1 );
      }
      else if ( kflag == 6 && (*relerr) <= relerr_save && abserr <= abserr_save )
      {
        cerr << "\n";
        cerr << "R8_GoldSaxRKF45ODESolve - Fatal error!\n";
        cerr << "  KFLAG = 6 and\n";
        cerr << "  RELERR <= RELERR_SAVE and\n";
        cerr << "  ABSERR <= ABSERR_SAVE\n";
        exit ( 1 );
      }
    }
//
//  FLAG = 3, 4, 5, 6, 7 or 8.
//
    else
    {
      if ( flag_return == 3 )
      {
        flag_return = flag_save;
        if ( kflag == 3 )
        {
          mflag = abs ( flag_return );
        }
      }
      else if ( flag_return == 4 )
      {
        nfe = 0;
        flag_return = flag_save;
        if ( kflag == 3 )
        {
          mflag = abs ( flag_return );
        }
      }
      else if ( flag_return == 5 && 0.0 < abserr )
      {
        flag_return = flag_save;
        if ( kflag == 3 )
        {
          mflag = abs ( flag_return );
        }
      }
//
//  Integration cannot be continued because the user did not respond to
//  the instructions pertaining to FLAG = 5, 6, 7 or 8.
//
      else
      {
        cerr << "\n";
        cerr << "R8_GoldSaxRKF45ODESolve - Fatal error!\n";
        cerr << "  Integration cannot be continued.\n";
        cerr << "  The user did not respond to the output\n";
        cerr << "  value FLAG = 5, 6, 7, or 8.\n";
        exit ( 1 );
      }
    }
  }
//
//  Save the input value of FLAG.  
//  Set the continuation flag KFLAG for subsequent input checking.
//
  flag_save = flag_return;
  kflag = 0;
//
//  Save RELERR and ABSERR for checking input on subsequent calls.
//
  relerr_save = (*relerr);
  abserr_save = abserr;
//
//  Restrict the relative error tolerance to be at least 
//
//    2*EPS+REMIN 
//
//  to avoid limiting precision difficulties arising from impossible 
//  accuracy requests.
//
  relerr_min = 2.0 * r8_epsilon ( ) + remin;
//
//  Is the relative error tolerance too small?
//
  if ( (*relerr) < relerr_min )
  {
    (*relerr) = relerr_min;
    kflag = 3;
    flag_return = 3;
    return flag_return;
  }

  dt = tout - *t;

  f1 = new double[neqn];
  f2 = new double[neqn];
  f3 = new double[neqn];
  f4 = new double[neqn];
  f5 = new double[neqn];

  if ( mflag == 1 )
  {
    init = 0;
    kop = 0;
    f ( *t, y, yp );
    nfe = 1;

    if ( *t == tout )
    {
      flag_return = 2;
      return flag_return;
    }

  }

  if ( init == 0 )
  {
    init = 1;
    h = r8_abs ( dt );
    toln = 0.0;

    for ( k = 0; k < neqn; k++ )
    {
      tol = (*relerr) * r8_abs ( y[k] ) + abserr;
      if ( 0.0 < tol )
      {
        toln = tol;
        ypk = r8_abs ( yp[k] );
        if ( tol < ypk * pow ( h, 5 ) )
        {
          h = pow ( ( tol / ypk ), 0.2 );
        }
      }
    }

    if ( toln <= 0.0 )
    {
      h = 0.0;
    }

    h = r8_max ( h, 26.0 * eps * r8_max ( r8_abs ( *t ), r8_abs ( dt ) ) );

    if ( flag_return < 0 )
    {
      flag_save = -2;
    }
    else
    {
      flag_save = 2;
    }
  }
//
//  Set stepsize for integration in the direction from T to TOUT.
//
  h = r8_sign ( dt ) * r8_abs ( h );
//
//  Test to see if too may output points are being requested.
//
  if ( 2.0 * r8_abs ( dt ) <= r8_abs ( h ) )
  {
    kop = kop + 1;
  }
//
//  Unnecessary frequency of output.
//
  if ( kop == 100 )
  {
    kop = 0;
    delete [] f1;
    delete [] f2;
    delete [] f3;
    delete [] f4;
    delete [] f5;
    flag_return = 7;
    return flag_return;
  }
//
//  If we are too close to the output point, then simply extrapolate and return.
//
  if ( r8_abs ( dt ) <= 26.0 * eps * r8_abs ( *t ) )
  {
    *t = tout;
    for ( i = 0; i < neqn; i++ )
    {
      y[i] = y[i] + dt * yp[i];
    }
    f ( *t, y, yp );
    nfe = nfe + 1;

    delete [] f1;
    delete [] f2;
    delete [] f3;
    delete [] f4;
    delete [] f5;
    flag_return = 2;
    return flag_return;
  }

  output = false;

  scale = 2.0 / (*relerr);
  ae = scale * abserr;

  for ( ; ; )
  {
    hfaild = false;

    hmin = 26.0 * eps * r8_abs ( *t );

    dt = tout - *t;

    if ( 2.0 * r8_abs ( h ) <= r8_abs ( dt ) )
    {
    }
    else

    {
      if ( r8_abs ( dt ) <= r8_abs ( h ) )
      {
        output = true;
        h = dt;
      }
      else
      {
        h = 0.5 * dt;
      }

    }

    for ( ; ; )
    {

      if ( MAXNFE < nfe )
      {
        kflag = 4;
        delete [] f1;
        delete [] f2;
        delete [] f3;
        delete [] f4;
        delete [] f5;
        flag_return = 4;
        return flag_return;
      }

      r8_fehl ( f, neqn, y, *t, h, yp, f1, f2, f3, f4, f5, f1 );
      nfe = nfe + 5;

      eeoet = 0.0;
 
      for ( k = 0; k < neqn; k++ )
      {
        et = r8_abs ( y[k] ) + r8_abs ( f1[k] ) + ae;

        if ( et <= 0.0 )
        {
          delete [] f1;
          delete [] f2;
          delete [] f3;
          delete [] f4;
          delete [] f5;
          flag_return = 5;
          return flag_return;
        }

        ee = r8_abs 
        ( ( -2090.0 * yp[k] 
          + ( 21970.0 * f3[k] - 15048.0 * f4[k] ) 
          ) 
        + ( 22528.0 * f2[k] - 27360.0 * f5[k] ) 
        );

        eeoet = r8_max ( eeoet, ee / et );

      }

      esttol = r8_abs ( h ) * eeoet * scale / 752400.0;

      if ( esttol <= 1.0 )
      {
        break;
      }

      hfaild = true;
      output = false;

      if ( esttol < 59049.0 )
      {
        s = 0.9 / pow ( esttol, 0.2 );
      }
      else
      {
        s = 0.1;
      }

      h = s * h;

      if ( r8_abs ( h ) < hmin )
      {
        kflag = 6;
        delete [] f1;
        delete [] f2;
        delete [] f3;
        delete [] f4;
        delete [] f5;
        flag_return = 6;
        return flag_return;
      }

    }

    *t = *t + h;
    for ( i = 0; i < neqn; i++ )
    {
      y[i] = f1[i];
    }
    f ( *t, y, yp );
    nfe = nfe + 1;

    if ( 0.0001889568 < esttol )
    {
      s = 0.9 / pow ( esttol, 0.2 );
    }
    else
    {
      s = 5.0;
    }

    if ( hfaild )
    {
      s = r8_min ( s, 1.0 );
    }

    h = r8_sign ( h ) * r8_max ( s * r8_abs ( h ), hmin );

    if ( output )
    {
      *t = tout;
      delete [] f1;
      delete [] f2;
      delete [] f3;
      delete [] f4;
      delete [] f5;
      flag_return = 2;
      return flag_return;
    }

    if ( flag_return <= 0 )
    {
      delete [] f1;
      delete [] f2;
      delete [] f3;
      delete [] f4;
      delete [] f5;
      flag_return = -2;
      return flag_return;
    }

  }
# undef MAXNFE
}

double r8_sign ( double x )

{
  if ( x < 0.0 )
  {
    return ( -1.0 );
  } 
  else
  {
    return ( +1.0 );
  }
}

void timestamp ( )

{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
