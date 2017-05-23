//-------------------------------------------------------------------------//
//                                                                         //
//  This benchmark is a serial C version of the NPB SP code. This C        //
//  version is developed by the Center for Manycore Programming at Seoul   //
//  National University and derived from the serial Fortran versions in    //
//  "NPB3.3-SER" developed by NAS.                                         //
//                                                                         //
//  Permission to use, copy, distribute and modify this software for any   //
//  purpose with or without fee is hereby granted. This software is        //
//  provided "as is" without express or implied warranty.                  //
//                                                                         //
//  Information on NPB 3.3, including the technical report, the original   //
//  specifications, source code, results and information on how to submit  //
//  new results, is available at:                                          //
//                                                                         //
//           http://www.nas.nasa.gov/Software/NPB/                         //
//                                                                         //
//  Send comments or suggestions for this C version to cmp@aces.snu.ac.kr  //
//                                                                         //
//          Center for Manycore Programming                                //
//          School of Computer Science and Engineering                     //
//          Seoul National University                                      //
//          Seoul 151-744, Korea                                           //
//                                                                         //
//          E-mail:  cmp@aces.snu.ac.kr                                    //
//                                                                         //
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
// Authors: Sangmin Seo, Jungwon Kim, Jun Lee, Jeongho Nah, Gangwon Jo,    //
//          and Jaejin Lee                                                 //
//-------------------------------------------------------------------------//

#include "header.h"

void adi()
{
  compute_rhs();

  txinvr();

  x_solve();

  y_solve();

  z_solve();

  add();
}


void ninvr()
{
  int i, j, k;
  double r1, r2, r3, r4, r5, t1, t2;

  if (timeron) timer_start(t_ninvr);
  for (k = 1; k <= nz2; k++) {
    for (j = 1; j <= ny2; j++) {
      for (i = 1; i <= nx2; i++) {
        r1 = rhs[k][j][i][0];
        r2 = rhs[k][j][i][1];
        r3 = rhs[k][j][i][2];
        r4 = rhs[k][j][i][3];
        r5 = rhs[k][j][i][4];

        t1 = bt * r3;
        t2 = 0.5 * ( r4 + r5 );

        rhs[k][j][i][0] = -r2;
        rhs[k][j][i][1] =  r1;
        rhs[k][j][i][2] = bt * ( r4 - r5 );
        rhs[k][j][i][3] = -t1 + t2;
        rhs[k][j][i][4] =  t1 + t2;
      }
    }
  }
  if (timeron) timer_stop(t_ninvr);
}

void pinvr()
{
  int i, j, k;
  double r1, r2, r3, r4, r5, t1, t2;

  if (timeron) timer_start(t_pinvr);
  for (k = 1; k <= nz2; k++) {
    for (j = 1; j <= ny2; j++) {
      for (i = 1; i <= nx2; i++) {
        r1 = rhs[k][j][i][0];
        r2 = rhs[k][j][i][1];
        r3 = rhs[k][j][i][2];
        r4 = rhs[k][j][i][3];
        r5 = rhs[k][j][i][4];

        t1 = bt * r1;
        t2 = 0.5 * ( r4 + r5 );

        rhs[k][j][i][0] =  bt * ( r4 - r5 );
        rhs[k][j][i][1] = -r3;
        rhs[k][j][i][2] =  r2;
        rhs[k][j][i][3] = -t1 + t2;
        rhs[k][j][i][4] =  t1 + t2;
      }
    }
  }
  if (timeron) timer_stop(t_pinvr);
}

void tzetar()
{
  int i, j, k;
  double t1, t2, t3, ac, xvel, yvel, zvel, r1, r2, r3, r4, r5;
  double btuz, ac2u, uzik1;

  if (timeron) timer_start(t_tzetar);
  for (k = 1; k <= nz2; k++) {
    for (j = 1; j <= ny2; j++) {
      for (i = 1; i <= nx2; i++) {
        xvel = us[k][j][i];
        yvel = vs[k][j][i];
        zvel = ws[k][j][i];
        ac   = speed[k][j][i];

        ac2u = ac*ac;

        r1 = rhs[k][j][i][0];
        r2 = rhs[k][j][i][1];
        r3 = rhs[k][j][i][2];
        r4 = rhs[k][j][i][3];
        r5 = rhs[k][j][i][4];

        uzik1 = u[k][j][i][0];
        btuz  = bt * uzik1;

        t1 = btuz/ac * (r4 + r5);
        t2 = r3 + t1;
        t3 = btuz * (r4 - r5);

        rhs[k][j][i][0] = t2;
        rhs[k][j][i][1] = -uzik1*r2 + xvel*t2;
        rhs[k][j][i][2] =  uzik1*r1 + yvel*t2;
        rhs[k][j][i][3] =  zvel*t2  + t3;
        rhs[k][j][i][4] =  uzik1*(-xvel*r2 + yvel*r1) +
                           qs[k][j][i]*t2 + c2iv*ac2u*t1 + zvel*t3;
      }
    }
  }
  if (timeron) timer_stop(t_tzetar);
}

void x_solve()
{
  int i, j, k, i1, i2, m;
  double ru1, fac1, fac2;

  if (timeron) timer_start(t_xsolve);
  for (k = 1; k <= nz2; k++) {
    lhsinit(nx2+1, ny2);

    //---------------------------------------------------------------------
    // Computes the left hand side for the three x-factors
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    // first fill the lhs for the u-eigenvalue
    //---------------------------------------------------------------------
    for (j = 1; j <= ny2; j++) {
      for (i = 0; i <= grid_points[0]-1; i++) {
        ru1 = c3c4*rho_i[k][j][i];
        cv[i] = us[k][j][i];
        rhon[i] = max(max(dx2+con43*ru1,dx5+c1c5*ru1), max(dxmax+ru1,dx1));
      }

      for (i = 1; i <= nx2; i++) {
        lhs[j][i][0] =  0.0;
        lhs[j][i][1] = -dttx2 * cv[i-1] - dttx1 * rhon[i-1];
        lhs[j][i][2] =  1.0 + c2dttx1 * rhon[i];
        lhs[j][i][3] =  dttx2 * cv[i+1] - dttx1 * rhon[i+1];
        lhs[j][i][4] =  0.0;
      }
    }

    //---------------------------------------------------------------------
    // add fourth order dissipation
    //---------------------------------------------------------------------
    for (j = 1; j <= ny2; j++) {
      i = 1;
      lhs[j][i][2] = lhs[j][i][2] + comz5;
      lhs[j][i][3] = lhs[j][i][3] - comz4;
      lhs[j][i][4] = lhs[j][i][4] + comz1;

      lhs[j][i+1][1] = lhs[j][i+1][1] - comz4;
      lhs[j][i+1][2] = lhs[j][i+1][2] + comz6;
      lhs[j][i+1][3] = lhs[j][i+1][3] - comz4;
      lhs[j][i+1][4] = lhs[j][i+1][4] + comz1;
    }

    for (j = 1; j <= ny2; j++) {
      for (i = 3; i <= grid_points[0]-4; i++) {
        lhs[j][i][0] = lhs[j][i][0] + comz1;
        lhs[j][i][1] = lhs[j][i][1] - comz4;
        lhs[j][i][2] = lhs[j][i][2] + comz6;
        lhs[j][i][3] = lhs[j][i][3] - comz4;
        lhs[j][i][4] = lhs[j][i][4] + comz1;
      }
    }

    for (j = 1; j <= ny2; j++) {
      i = grid_points[0]-3;
      lhs[j][i][0] = lhs[j][i][0] + comz1;
      lhs[j][i][1] = lhs[j][i][1] - comz4;
      lhs[j][i][2] = lhs[j][i][2] + comz6;
      lhs[j][i][3] = lhs[j][i][3] - comz4;

      lhs[j][i+1][0] = lhs[j][i+1][0] + comz1;
      lhs[j][i+1][1] = lhs[j][i+1][1] - comz4;
      lhs[j][i+1][2] = lhs[j][i+1][2] + comz5;
    }

    //---------------------------------------------------------------------
    // subsequently, fill the other factors (u+c), (u-c) by adding to
    // the first
    //---------------------------------------------------------------------
    for (j = 1; j <= ny2; j++) {
      for (i = 1; i <= nx2; i++) {
        lhsp[j][i][0] = lhs[j][i][0];
        lhsp[j][i][1] = lhs[j][i][1] - dttx2 * speed[k][j][i-1];
        lhsp[j][i][2] = lhs[j][i][2];
        lhsp[j][i][3] = lhs[j][i][3] + dttx2 * speed[k][j][i+1];
        lhsp[j][i][4] = lhs[j][i][4];
        lhsm[j][i][0] = lhs[j][i][0];
        lhsm[j][i][1] = lhs[j][i][1] + dttx2 * speed[k][j][i-1];
        lhsm[j][i][2] = lhs[j][i][2];
        lhsm[j][i][3] = lhs[j][i][3] - dttx2 * speed[k][j][i+1];
        lhsm[j][i][4] = lhs[j][i][4];
      }
    }

    //---------------------------------------------------------------------
    // FORWARD ELIMINATION
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    // perform the Thomas algorithm; first, FORWARD ELIMINATION
    //---------------------------------------------------------------------
    for (j = 1; j <= ny2; j++) {
      for (i = 0; i <= grid_points[0]-3; i++) {
        i1 = i + 1;
        i2 = i + 2;
        fac1 = 1.0/lhs[j][i][2];
        lhs[j][i][3] = fac1*lhs[j][i][3];
        lhs[j][i][4] = fac1*lhs[j][i][4];
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
        }
        lhs[j][i1][2] = lhs[j][i1][2] - lhs[j][i1][1]*lhs[j][i][3];
        lhs[j][i1][3] = lhs[j][i1][3] - lhs[j][i1][1]*lhs[j][i][4];
        for (m = 0; m < 3; m++) {
          rhs[k][j][i1][m] = rhs[k][j][i1][m] - lhs[j][i1][1]*rhs[k][j][i][m];
        }
        lhs[j][i2][1] = lhs[j][i2][1] - lhs[j][i2][0]*lhs[j][i][3];
        lhs[j][i2][2] = lhs[j][i2][2] - lhs[j][i2][0]*lhs[j][i][4];
        for (m = 0; m < 3; m++) {
          rhs[k][j][i2][m] = rhs[k][j][i2][m] - lhs[j][i2][0]*rhs[k][j][i][m];
        }
      }
    }

    //---------------------------------------------------------------------
    // The last two rows in this grid block are a bit different,
    // since they for (not have two more rows available for the
    // elimination of off-diagonal entries
    //---------------------------------------------------------------------
    for (j = 1; j <= ny2; j++) {
      i  = grid_points[0]-2;
      i1 = grid_points[0]-1;
      fac1 = 1.0/lhs[j][i][2];
      lhs[j][i][3] = fac1*lhs[j][i][3];
      lhs[j][i][4] = fac1*lhs[j][i][4];
      for (m = 0; m < 3; m++) {
        rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
      }
      lhs[j][i1][2] = lhs[j][i1][2] - lhs[j][i1][1]*lhs[j][i][3];
      lhs[j][i1][3] = lhs[j][i1][3] - lhs[j][i1][1]*lhs[j][i][4];
      for (m = 0; m < 3; m++) {
        rhs[k][j][i1][m] = rhs[k][j][i1][m] - lhs[j][i1][1]*rhs[k][j][i][m];
      }

      //---------------------------------------------------------------------
      // scale the last row immediately
      //---------------------------------------------------------------------
      fac2 = 1.0/lhs[j][i1][2];
      for (m = 0; m < 3; m++) {
        rhs[k][j][i1][m] = fac2*rhs[k][j][i1][m];
      }
    }

    //---------------------------------------------------------------------
    // for (the u+c and the u-c factors
    //---------------------------------------------------------------------
    for (j = 1; j <= ny2; j++) {
      for (i = 0; i <= grid_points[0]-3; i++) {
        i1 = i + 1;
        i2 = i + 2;

        m = 3;
        fac1 = 1.0/lhsp[j][i][2];
        lhsp[j][i][3]    = fac1*lhsp[j][i][3];
        lhsp[j][i][4]    = fac1*lhsp[j][i][4];
        rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
        lhsp[j][i1][2]   = lhsp[j][i1][2] - lhsp[j][i1][1]*lhsp[j][i][3];
        lhsp[j][i1][3]   = lhsp[j][i1][3] - lhsp[j][i1][1]*lhsp[j][i][4];
        rhs[k][j][i1][m] = rhs[k][j][i1][m] - lhsp[j][i1][1]*rhs[k][j][i][m];
        lhsp[j][i2][1]   = lhsp[j][i2][1] - lhsp[j][i2][0]*lhsp[j][i][3];
        lhsp[j][i2][2]   = lhsp[j][i2][2] - lhsp[j][i2][0]*lhsp[j][i][4];
        rhs[k][j][i2][m] = rhs[k][j][i2][m] - lhsp[j][i2][0]*rhs[k][j][i][m];

        m = 4;
        fac1 = 1.0/lhsm[j][i][2];
        lhsm[j][i][3]    = fac1*lhsm[j][i][3];
        lhsm[j][i][4]    = fac1*lhsm[j][i][4];
        rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
        lhsm[j][i1][2]   = lhsm[j][i1][2] - lhsm[j][i1][1]*lhsm[j][i][3];
        lhsm[j][i1][3]   = lhsm[j][i1][3] - lhsm[j][i1][1]*lhsm[j][i][4];
        rhs[k][j][i1][m] = rhs[k][j][i1][m] - lhsm[j][i1][1]*rhs[k][j][i][m];
        lhsm[j][i2][1]   = lhsm[j][i2][1] - lhsm[j][i2][0]*lhsm[j][i][3];
        lhsm[j][i2][2]   = lhsm[j][i2][2] - lhsm[j][i2][0]*lhsm[j][i][4];
        rhs[k][j][i2][m] = rhs[k][j][i2][m] - lhsm[j][i2][0]*rhs[k][j][i][m];
      }
    }

    //---------------------------------------------------------------------
    // And again the last two rows separately
    //---------------------------------------------------------------------
    for (j = 1; j <= ny2; j++) {
      i  = grid_points[0]-2;
      i1 = grid_points[0]-1;

      m = 3;
      fac1 = 1.0/lhsp[j][i][2];
      lhsp[j][i][3]    = fac1*lhsp[j][i][3];
      lhsp[j][i][4]    = fac1*lhsp[j][i][4];
      rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
      lhsp[j][i1][2]   = lhsp[j][i1][2] - lhsp[j][i1][1]*lhsp[j][i][3];
      lhsp[j][i1][3]   = lhsp[j][i1][3] - lhsp[j][i1][1]*lhsp[j][i][4];
      rhs[k][j][i1][m] = rhs[k][j][i1][m] - lhsp[j][i1][1]*rhs[k][j][i][m];

      m = 4;
      fac1 = 1.0/lhsm[j][i][2];
      lhsm[j][i][3]    = fac1*lhsm[j][i][3];
      lhsm[j][i][4]    = fac1*lhsm[j][i][4];
      rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
      lhsm[j][i1][2]   = lhsm[j][i1][2] - lhsm[j][i1][1]*lhsm[j][i][3];
      lhsm[j][i1][3]   = lhsm[j][i1][3] - lhsm[j][i1][1]*lhsm[j][i][4];
      rhs[k][j][i1][m] = rhs[k][j][i1][m] - lhsm[j][i1][1]*rhs[k][j][i][m];

      //---------------------------------------------------------------------
      // Scale the last row immediately
      //---------------------------------------------------------------------
      rhs[k][j][i1][3] = rhs[k][j][i1][3]/lhsp[j][i1][2];
      rhs[k][j][i1][4] = rhs[k][j][i1][4]/lhsm[j][i1][2];
    }

    //---------------------------------------------------------------------
    // BACKSUBSTITUTION
    //---------------------------------------------------------------------
    for (j = 1; j <= ny2; j++) {
      i  = grid_points[0]-2;
      i1 = grid_points[0]-1;
      for (m = 0; m < 3; m++) {
        rhs[k][j][i][m] = rhs[k][j][i][m] - lhs[j][i][3]*rhs[k][j][i1][m];
      }

      rhs[k][j][i][3] = rhs[k][j][i][3] - lhsp[j][i][3]*rhs[k][j][i1][3];
      rhs[k][j][i][4] = rhs[k][j][i][4] - lhsm[j][i][3]*rhs[k][j][i1][4];
    }

    //---------------------------------------------------------------------
    // The first three factors
    //---------------------------------------------------------------------
    for (j = 1; j <= ny2; j++) {
      for (i = grid_points[0]-3; i >= 0; i--) {
        i1 = i + 1;
        i2 = i + 2;
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] = rhs[k][j][i][m] -
                            lhs[j][i][3]*rhs[k][j][i1][m] -
                            lhs[j][i][4]*rhs[k][j][i2][m];
        }

        //-------------------------------------------------------------------
        // And the remaining two
        //-------------------------------------------------------------------
        rhs[k][j][i][3] = rhs[k][j][i][3] -
                          lhsp[j][i][3]*rhs[k][j][i1][3] -
                          lhsp[j][i][4]*rhs[k][j][i2][3];
        rhs[k][j][i][4] = rhs[k][j][i][4] -
                          lhsm[j][i][3]*rhs[k][j][i1][4] -
                          lhsm[j][i][4]*rhs[k][j][i2][4];
      }
    }
  }

  //---------------------------------------------------------------------
  // Do the block-diagonal inversion
  //---------------------------------------------------------------------
  ninvr();
}

void y_solve()
{
  int i, j, k, j1, j2, m;
  double ru1, fac1, fac2;

  if (timeron) timer_start(t_ysolve);
  for (k = 1; k <= grid_points[2]-2; k++) {
    lhsinitj(ny2+1, nx2);

    //---------------------------------------------------------------------
    // Computes the left hand side for the three y-factors
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    // first fill the lhs for the u-eigenvalue
    //---------------------------------------------------------------------
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 0; j <= grid_points[1]-1; j++) {
        ru1 = c3c4*rho_i[k][j][i];
        cv[j] = vs[k][j][i];
        rhoq[j] = max(max(dy3+con43*ru1, dy5+c1c5*ru1), max(dymax+ru1, dy1));
      }

      for (j = 1; j <= grid_points[1]-2; j++) {
        lhs[j][i][0] =  0.0;
        lhs[j][i][1] = -dtty2 * cv[j-1] - dtty1 * rhoq[j-1];
        lhs[j][i][2] =  1.0 + c2dtty1 * rhoq[j];
        lhs[j][i][3] =  dtty2 * cv[j+1] - dtty1 * rhoq[j+1];
        lhs[j][i][4] =  0.0;
      }
    }

    //---------------------------------------------------------------------
    // add fourth order dissipation
    //---------------------------------------------------------------------
    for (i = 1; i <= grid_points[0]-2; i++) {
      j = 1;
      lhs[j][i][2] = lhs[j][i][2] + comz5;
      lhs[j][i][3] = lhs[j][i][3] - comz4;
      lhs[j][i][4] = lhs[j][i][4] + comz1;

      lhs[j+1][i][1] = lhs[j+1][i][1] - comz4;
      lhs[j+1][i][2] = lhs[j+1][i][2] + comz6;
      lhs[j+1][i][3] = lhs[j+1][i][3] - comz4;
      lhs[j+1][i][4] = lhs[j+1][i][4] + comz1;
    }

    for (j = 3; j <= grid_points[1]-4; j++) {
      for (i = 1; i <= grid_points[0]-2; i++) {
        lhs[j][i][0] = lhs[j][i][0] + comz1;
        lhs[j][i][1] = lhs[j][i][1] - comz4;
        lhs[j][i][2] = lhs[j][i][2] + comz6;
        lhs[j][i][3] = lhs[j][i][3] - comz4;
        lhs[j][i][4] = lhs[j][i][4] + comz1;
      }
    }

    for (i = 1; i <= grid_points[0]-2; i++) {
      j = grid_points[1]-3;
      lhs[j][i][0] = lhs[j][i][0] + comz1;
      lhs[j][i][1] = lhs[j][i][1] - comz4;
      lhs[j][i][2] = lhs[j][i][2] + comz6;
      lhs[j][i][3] = lhs[j][i][3] - comz4;

      lhs[j+1][i][0] = lhs[j+1][i][0] + comz1;
      lhs[j+1][i][1] = lhs[j+1][i][1] - comz4;
      lhs[j+1][i][2] = lhs[j+1][i][2] + comz5;
    }

    //---------------------------------------------------------------------
    // subsequently, for (the other two factors
    //---------------------------------------------------------------------
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (i = 1; i <= grid_points[0]-2; i++) {
        lhsp[j][i][0] = lhs[j][i][0];
        lhsp[j][i][1] = lhs[j][i][1] - dtty2 * speed[k][j-1][i];
        lhsp[j][i][2] = lhs[j][i][2];
        lhsp[j][i][3] = lhs[j][i][3] + dtty2 * speed[k][j+1][i];
        lhsp[j][i][4] = lhs[j][i][4];
        lhsm[j][i][0] = lhs[j][i][0];
        lhsm[j][i][1] = lhs[j][i][1] + dtty2 * speed[k][j-1][i];
        lhsm[j][i][2] = lhs[j][i][2];
        lhsm[j][i][3] = lhs[j][i][3] - dtty2 * speed[k][j+1][i];
        lhsm[j][i][4] = lhs[j][i][4];
      }
    }


    //---------------------------------------------------------------------
    // FORWARD ELIMINATION
    //---------------------------------------------------------------------
    for (j = 0; j <= grid_points[1]-3; j++) {
      j1 = j + 1;
      j2 = j + 2;
      for (i = 1; i <= grid_points[0]-2; i++) {
        fac1 = 1.0/lhs[j][i][2];
        lhs[j][i][3] = fac1*lhs[j][i][3];
        lhs[j][i][4] = fac1*lhs[j][i][4];
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
        }
        lhs[j1][i][2] = lhs[j1][i][2] - lhs[j1][i][1]*lhs[j][i][3];
        lhs[j1][i][3] = lhs[j1][i][3] - lhs[j1][i][1]*lhs[j][i][4];
        for (m = 0; m < 3; m++) {
          rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhs[j1][i][1]*rhs[k][j][i][m];
        }
        lhs[j2][i][1] = lhs[j2][i][1] - lhs[j2][i][0]*lhs[j][i][3];
        lhs[j2][i][2] = lhs[j2][i][2] - lhs[j2][i][0]*lhs[j][i][4];
        for (m = 0; m < 3; m++) {
          rhs[k][j2][i][m] = rhs[k][j2][i][m] - lhs[j2][i][0]*rhs[k][j][i][m];
        }
      }
    }

    //---------------------------------------------------------------------
    // The last two rows in this grid block are a bit different,
    // since they for (not have two more rows available for the
    // elimination of off-diagonal entries
    //---------------------------------------------------------------------
    j  = grid_points[1]-2;
    j1 = grid_points[1]-1;
    for (i = 1; i <= grid_points[0]-2; i++) {
      fac1 = 1.0/lhs[j][i][2];
      lhs[j][i][3] = fac1*lhs[j][i][3];
      lhs[j][i][4] = fac1*lhs[j][i][4];
      for (m = 0; m < 3; m++) {
        rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
      }
      lhs[j1][i][2] = lhs[j1][i][2] - lhs[j1][i][1]*lhs[j][i][3];
      lhs[j1][i][3] = lhs[j1][i][3] - lhs[j1][i][1]*lhs[j][i][4];
      for (m = 0; m < 3; m++) {
        rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhs[j1][i][1]*rhs[k][j][i][m];
      }
      //---------------------------------------------------------------------
      // scale the last row immediately
      //---------------------------------------------------------------------
      fac2 = 1.0/lhs[j1][i][2];
      for (m = 0; m < 3; m++) {
        rhs[k][j1][i][m] = fac2*rhs[k][j1][i][m];
      }
    }

    //---------------------------------------------------------------------
    // for (the u+c and the u-c factors
    //---------------------------------------------------------------------
    for (j = 0; j <= grid_points[1]-3; j++) {
      j1 = j + 1;
      j2 = j + 2;
      for (i = 1; i <= grid_points[0]-2; i++) {
        m = 3;
        fac1 = 1.0/lhsp[j][i][2];
        lhsp[j][i][3]    = fac1*lhsp[j][i][3];
        lhsp[j][i][4]    = fac1*lhsp[j][i][4];
        rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
        lhsp[j1][i][2]   = lhsp[j1][i][2] - lhsp[j1][i][1]*lhsp[j][i][3];
        lhsp[j1][i][3]   = lhsp[j1][i][3] - lhsp[j1][i][1]*lhsp[j][i][4];
        rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhsp[j1][i][1]*rhs[k][j][i][m];
        lhsp[j2][i][1]   = lhsp[j2][i][1] - lhsp[j2][i][0]*lhsp[j][i][3];
        lhsp[j2][i][2]   = lhsp[j2][i][2] - lhsp[j2][i][0]*lhsp[j][i][4];
        rhs[k][j2][i][m] = rhs[k][j2][i][m] - lhsp[j2][i][0]*rhs[k][j][i][m];

        m = 4;
        fac1 = 1.0/lhsm[j][i][2];
        lhsm[j][i][3]    = fac1*lhsm[j][i][3];
        lhsm[j][i][4]    = fac1*lhsm[j][i][4];
        rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
        lhsm[j1][i][2]   = lhsm[j1][i][2] - lhsm[j1][i][1]*lhsm[j][i][3];
        lhsm[j1][i][3]   = lhsm[j1][i][3] - lhsm[j1][i][1]*lhsm[j][i][4];
        rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhsm[j1][i][1]*rhs[k][j][i][m];
        lhsm[j2][i][1]   = lhsm[j2][i][1] - lhsm[j2][i][0]*lhsm[j][i][3];
        lhsm[j2][i][2]   = lhsm[j2][i][2] - lhsm[j2][i][0]*lhsm[j][i][4];
        rhs[k][j2][i][m] = rhs[k][j2][i][m] - lhsm[j2][i][0]*rhs[k][j][i][m];
      }
    }

    //---------------------------------------------------------------------
    // And again the last two rows separately
    //---------------------------------------------------------------------
    j  = grid_points[1]-2;
    j1 = grid_points[1]-1;
    for (i = 1; i <= grid_points[0]-2; i++) {
      m = 3;
      fac1 = 1.0/lhsp[j][i][2];
      lhsp[j][i][3]    = fac1*lhsp[j][i][3];
      lhsp[j][i][4]    = fac1*lhsp[j][i][4];
      rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
      lhsp[j1][i][2]   = lhsp[j1][i][2] - lhsp[j1][i][1]*lhsp[j][i][3];
      lhsp[j1][i][3]   = lhsp[j1][i][3] - lhsp[j1][i][1]*lhsp[j][i][4];
      rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhsp[j1][i][1]*rhs[k][j][i][m];

      m = 4;
      fac1 = 1.0/lhsm[j][i][2];
      lhsm[j][i][3]    = fac1*lhsm[j][i][3];
      lhsm[j][i][4]    = fac1*lhsm[j][i][4];
      rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
      lhsm[j1][i][2]   = lhsm[j1][i][2] - lhsm[j1][i][1]*lhsm[j][i][3];
      lhsm[j1][i][3]   = lhsm[j1][i][3] - lhsm[j1][i][1]*lhsm[j][i][4];
      rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhsm[j1][i][1]*rhs[k][j][i][m];

      //---------------------------------------------------------------------
      // Scale the last row immediately
      //---------------------------------------------------------------------
      rhs[k][j1][i][3]   = rhs[k][j1][i][3]/lhsp[j1][i][2];
      rhs[k][j1][i][4]   = rhs[k][j1][i][4]/lhsm[j1][i][2];
    }


    //---------------------------------------------------------------------
    // BACKSUBSTITUTION
    //---------------------------------------------------------------------
    j  = grid_points[1]-2;
    j1 = grid_points[1]-1;
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (m = 0; m < 3; m++) {
        rhs[k][j][i][m] = rhs[k][j][i][m] - lhs[j][i][3]*rhs[k][j1][i][m];
      }

      rhs[k][j][i][3] = rhs[k][j][i][3] - lhsp[j][i][3]*rhs[k][j1][i][3];
      rhs[k][j][i][4] = rhs[k][j][i][4] - lhsm[j][i][3]*rhs[k][j1][i][4];
    }

    //---------------------------------------------------------------------
    // The first three factors
    //---------------------------------------------------------------------
    for (j = grid_points[1]-3; j >= 0; j--) {
      j1 = j + 1;
      j2 = j + 2;
      for (i = 1; i <= grid_points[0]-2; i++) {
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] = rhs[k][j][i][m] -
                            lhs[j][i][3]*rhs[k][j1][i][m] -
                            lhs[j][i][4]*rhs[k][j2][i][m];
        }

        //-------------------------------------------------------------------
        // And the remaining two
        //-------------------------------------------------------------------
        rhs[k][j][i][3] = rhs[k][j][i][3] -
                          lhsp[j][i][3]*rhs[k][j1][i][3] -
                          lhsp[j][i][4]*rhs[k][j2][i][3];
        rhs[k][j][i][4] = rhs[k][j][i][4] -
                          lhsm[j][i][3]*rhs[k][j1][i][4] -
                          lhsm[j][i][4]*rhs[k][j2][i][4];
      }
    }
  }

  pinvr();
}

void z_solve()
{
  int i, j, k, k1, k2, m;
  double ru1, fac1, fac2;

  for (j = 1; j <= ny2; j++) {
    lhsinitj(nz2+1, nx2);

    //---------------------------------------------------------------------
    // Computes the left hand side for the three z-factors
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    // first fill the lhs for the u-eigenvalue
    //---------------------------------------------------------------------
    for (i = 1; i <= nx2; i++) {
      for (k = 0; k <= nz2+1; k++) {
        ru1 = c3c4*rho_i[k][j][i];
        cv[k] = ws[k][j][i];
        rhos[k] = max(max(dz4+con43*ru1, dz5+c1c5*ru1), max(dzmax+ru1, dz1));
      }

      for (k = 1; k <= nz2; k++) {
        lhs[k][i][0] =  0.0;
        lhs[k][i][1] = -dttz2 * cv[k-1] - dttz1 * rhos[k-1];
        lhs[k][i][2] =  1.0 + c2dttz1 * rhos[k];
        lhs[k][i][3] =  dttz2 * cv[k+1] - dttz1 * rhos[k+1];
        lhs[k][i][4] =  0.0;
      }
    }

    //---------------------------------------------------------------------
    // add fourth order dissipation
    //---------------------------------------------------------------------
    for (i = 1; i <= nx2; i++) {
      k = 1;
      lhs[k][i][2] = lhs[k][i][2] + comz5;
      lhs[k][i][3] = lhs[k][i][3] - comz4;
      lhs[k][i][4] = lhs[k][i][4] + comz1;

      k = 2;
      lhs[k][i][1] = lhs[k][i][1] - comz4;
      lhs[k][i][2] = lhs[k][i][2] + comz6;
      lhs[k][i][3] = lhs[k][i][3] - comz4;
      lhs[k][i][4] = lhs[k][i][4] + comz1;
    }

    for (k = 3; k <= nz2-2; k++) {
      for (i = 1; i <= nx2; i++) {
        lhs[k][i][0] = lhs[k][i][0] + comz1;
        lhs[k][i][1] = lhs[k][i][1] - comz4;
        lhs[k][i][2] = lhs[k][i][2] + comz6;
        lhs[k][i][3] = lhs[k][i][3] - comz4;
        lhs[k][i][4] = lhs[k][i][4] + comz1;
      }
    }

    for (i = 1; i <= nx2; i++) {
      k = nz2-1;
      lhs[k][i][0] = lhs[k][i][0] + comz1;
      lhs[k][i][1] = lhs[k][i][1] - comz4;
      lhs[k][i][2] = lhs[k][i][2] + comz6;
      lhs[k][i][3] = lhs[k][i][3] - comz4;

      k = nz2;
      lhs[k][i][0] = lhs[k][i][0] + comz1;
      lhs[k][i][1] = lhs[k][i][1] - comz4;
      lhs[k][i][2] = lhs[k][i][2] + comz5;
    }

    //---------------------------------------------------------------------
    // subsequently, fill the other factors (u+c), (u-c)
    //---------------------------------------------------------------------
    for (k = 1; k <= nz2; k++) {
      for (i = 1; i <= nx2; i++) {
        lhsp[k][i][0] = lhs[k][i][0];
        lhsp[k][i][1] = lhs[k][i][1] - dttz2 * speed[k-1][j][i];
        lhsp[k][i][2] = lhs[k][i][2];
        lhsp[k][i][3] = lhs[k][i][3] + dttz2 * speed[k+1][j][i];
        lhsp[k][i][4] = lhs[k][i][4];
        lhsm[k][i][0] = lhs[k][i][0];
        lhsm[k][i][1] = lhs[k][i][1] + dttz2 * speed[k-1][j][i];
        lhsm[k][i][2] = lhs[k][i][2];
        lhsm[k][i][3] = lhs[k][i][3] - dttz2 * speed[k+1][j][i];
        lhsm[k][i][4] = lhs[k][i][4];
      }
    }


    //---------------------------------------------------------------------
    // FORWARD ELIMINATION
    //---------------------------------------------------------------------
    for (k = 0; k <= grid_points[2]-3; k++) {
      k1 = k + 1;
      k2 = k + 2;
      for (i = 1; i <= nx2; i++) {
        fac1 = 1.0/lhs[k][i][2];
        lhs[k][i][3] = fac1*lhs[k][i][3];
        lhs[k][i][4] = fac1*lhs[k][i][4];
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
        }
        lhs[k1][i][2] = lhs[k1][i][2] - lhs[k1][i][1]*lhs[k][i][3];
        lhs[k1][i][3] = lhs[k1][i][3] - lhs[k1][i][1]*lhs[k][i][4];
        for (m = 0; m < 3; m++) {
          rhs[k1][j][i][m] = rhs[k1][j][i][m] - lhs[k1][i][1]*rhs[k][j][i][m];
        }
        lhs[k2][i][1] = lhs[k2][i][1] - lhs[k2][i][0]*lhs[k][i][3];
        lhs[k2][i][2] = lhs[k2][i][2] - lhs[k2][i][0]*lhs[k][i][4];
        for (m = 0; m < 3; m++) {
          rhs[k2][j][i][m] = rhs[k2][j][i][m] - lhs[k2][i][0]*rhs[k][j][i][m];
        }
      }
    }

    //---------------------------------------------------------------------
    // The last two rows in this grid block are a bit different,
    // since they for (not have two more rows available for the
    // elimination of off-diagonal entries
    //---------------------------------------------------------------------
    k  = grid_points[2]-2;
    k1 = grid_points[2]-1;
    for (i = 1; i <= nx2; i++) {
      fac1 = 1.0/lhs[k][i][2];
      lhs[k][i][3] = fac1*lhs[k][i][3];
      lhs[k][i][4] = fac1*lhs[k][i][4];
      for (m = 0; m < 3; m++) {
        rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
      }
      lhs[k1][i][2] = lhs[k1][i][2] - lhs[k1][i][1]*lhs[k][i][3];
      lhs[k1][i][3] = lhs[k1][i][3] - lhs[k1][i][1]*lhs[k][i][4];
      for (m = 0; m < 3; m++) {
        rhs[k1][j][i][m] = rhs[k1][j][i][m] - lhs[k1][i][1]*rhs[k][j][i][m];
      }

      //---------------------------------------------------------------------
      // scale the last row immediately
      //---------------------------------------------------------------------
      fac2 = 1.0/lhs[k1][i][2];
      for (m = 0; m < 3; m++) {
        rhs[k1][j][i][m] = fac2*rhs[k1][j][i][m];
      }
    }

    //---------------------------------------------------------------------
    // for (the u+c and the u-c factors
    //---------------------------------------------------------------------
    for (k = 0; k <= grid_points[2]-3; k++) {
      k1 = k + 1;
      k2 = k + 2;
      for (i = 1; i <= nx2; i++) {
        m = 3;
        fac1 = 1.0/lhsp[k][i][2];
        lhsp[k][i][3]    = fac1*lhsp[k][i][3];
        lhsp[k][i][4]    = fac1*lhsp[k][i][4];
        rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
        lhsp[k1][i][2]   = lhsp[k1][i][2] - lhsp[k1][i][1]*lhsp[k][i][3];
        lhsp[k1][i][3]   = lhsp[k1][i][3] - lhsp[k1][i][1]*lhsp[k][i][4];
        rhs[k1][j][i][m] = rhs[k1][j][i][m] - lhsp[k1][i][1]*rhs[k][j][i][m];
        lhsp[k2][i][1]   = lhsp[k2][i][1] - lhsp[k2][i][0]*lhsp[k][i][3];
        lhsp[k2][i][2]   = lhsp[k2][i][2] - lhsp[k2][i][0]*lhsp[k][i][4];
        rhs[k2][j][i][m] = rhs[k2][j][i][m] - lhsp[k2][i][0]*rhs[k][j][i][m];

        m = 4;
        fac1 = 1.0/lhsm[k][i][2];
        lhsm[k][i][3]    = fac1*lhsm[k][i][3];
        lhsm[k][i][4]    = fac1*lhsm[k][i][4];
        rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
        lhsm[k1][i][2]   = lhsm[k1][i][2] - lhsm[k1][i][1]*lhsm[k][i][3];
        lhsm[k1][i][3]   = lhsm[k1][i][3] - lhsm[k1][i][1]*lhsm[k][i][4];
        rhs[k1][j][i][m] = rhs[k1][j][i][m] - lhsm[k1][i][1]*rhs[k][j][i][m];
        lhsm[k2][i][1]   = lhsm[k2][i][1] - lhsm[k2][i][0]*lhsm[k][i][3];
        lhsm[k2][i][2]   = lhsm[k2][i][2] - lhsm[k2][i][0]*lhsm[k][i][4];
        rhs[k2][j][i][m] = rhs[k2][j][i][m] - lhsm[k2][i][0]*rhs[k][j][i][m];
      }
    }

    //---------------------------------------------------------------------
    // And again the last two rows separately
    //---------------------------------------------------------------------
    k  = grid_points[2]-2;
    k1 = grid_points[2]-1;
    for (i = 1; i <= nx2; i++) {
      m = 3;
      fac1 = 1.0/lhsp[k][i][2];
      lhsp[k][i][3]    = fac1*lhsp[k][i][3];
      lhsp[k][i][4]    = fac1*lhsp[k][i][4];
      rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
      lhsp[k1][i][2]   = lhsp[k1][i][2] - lhsp[k1][i][1]*lhsp[k][i][3];
      lhsp[k1][i][3]   = lhsp[k1][i][3] - lhsp[k1][i][1]*lhsp[k][i][4];
      rhs[k1][j][i][m] = rhs[k1][j][i][m] - lhsp[k1][i][1]*rhs[k][j][i][m];

      m = 4;
      fac1 = 1.0/lhsm[k][i][2];
      lhsm[k][i][3]    = fac1*lhsm[k][i][3];
      lhsm[k][i][4]    = fac1*lhsm[k][i][4];
      rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
      lhsm[k1][i][2]   = lhsm[k1][i][2] - lhsm[k1][i][1]*lhsm[k][i][3];
      lhsm[k1][i][3]   = lhsm[k1][i][3] - lhsm[k1][i][1]*lhsm[k][i][4];
      rhs[k1][j][i][m] = rhs[k1][j][i][m] - lhsm[k1][i][1]*rhs[k][j][i][m];

      //---------------------------------------------------------------------
      // Scale the last row immediately (some of this is overkill
      // if this is the last cell)
      //---------------------------------------------------------------------
      rhs[k1][j][i][3] = rhs[k1][j][i][3]/lhsp[k1][i][2];
      rhs[k1][j][i][4] = rhs[k1][j][i][4]/lhsm[k1][i][2];
    }


    //---------------------------------------------------------------------
    // BACKSUBSTITUTION
    //---------------------------------------------------------------------
    k  = grid_points[2]-2;
    k1 = grid_points[2]-1;
    for (i = 1; i <= nx2; i++) {
      for (m = 0; m < 3; m++) {
        rhs[k][j][i][m] = rhs[k][j][i][m] - lhs[k][i][3]*rhs[k1][j][i][m];
      }

      rhs[k][j][i][3] = rhs[k][j][i][3] - lhsp[k][i][3]*rhs[k1][j][i][3];
      rhs[k][j][i][4] = rhs[k][j][i][4] - lhsm[k][i][3]*rhs[k1][j][i][4];
    }

    //---------------------------------------------------------------------
    // Whether or not this is the last processor, we always have
    // to complete the back-substitution
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    // The first three factors
    //---------------------------------------------------------------------
    for (k = grid_points[2]-3; k >= 0; k--) {
      k1 = k + 1;
      k2 = k + 2;
      for (i = 1; i <= nx2; i++) {
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] = rhs[k][j][i][m] -
                            lhs[k][i][3]*rhs[k1][j][i][m] -
                            lhs[k][i][4]*rhs[k2][j][i][m];
        }

        //-------------------------------------------------------------------
        // And the remaining two
        //-------------------------------------------------------------------
        rhs[k][j][i][3] = rhs[k][j][i][3] -
                          lhsp[k][i][3]*rhs[k1][j][i][3] -
                          lhsp[k][i][4]*rhs[k2][j][i][3];
        rhs[k][j][i][4] = rhs[k][j][i][4] -
                          lhsm[k][i][3]*rhs[k1][j][i][4] -
                          lhsm[k][i][4]*rhs[k2][j][i][4];
      }
    }
  }

  tzetar();
}
