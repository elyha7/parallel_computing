#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "mycom.h"
#include "mynet.h"
#include "myio.h"
#include "myprog.h"
int np , mp , nl , ier , lp ;
int mp_l , mp_r ;
char pname [MPI_MAX_PROCESSOR_NAME];
char vname [10] = "zadacha" ;
char sname [20];
MPI_Status status ;
union_t buf ;
double tick , t1 , t2 , t3 ;
double alpha = 0.1;
double omega = 20;
const double MAX_K = 7.0;
FILE * Fi = NULL ;
FILE * Fo = NULL ;
int nx , ntp , ntm , ntv , nyv ;
double xa , xb , xk , r0 , q0 , u0 , u1 ;
double theta , tau1 , tmax , epst ;
double tv , omg1 , gt ;
double f(double x,double t)
{
    double ret;
    ret=2*exp(-(x-0.5)*(x-0.5)/(alpha*alpha))*(1-exp(omega*t));
    return ret;
    return 0;
}
double k ( double x , double u ) ;
double k ( double x , double u ) {
    double u2 = u * u ;
    return 1+5* u2 ;
}
double k1 ( double u ) ;
double k1 ( double u ) {
    return 10* u ;
}
double q ( double u ) ;
double q ( double u ) {
    double su = cos ( u ) ;
    return su * su ;
}
double q1 ( double u ) ;
double q1 ( double u ) {
    return -sin (2* u ) ;
}
double g0 ( double x ) ;
double g0 ( double x ) {
    return 0 ;
}
double g1 ( double t ) ;
double g1 ( double t ) {
    return 0 ;
}
double g2 ( double t ) ;
double g2 ( double t ) {
    return 0 ;
}
int main ( int argc , char * argv [])
{
    int i , j , ii , i1 , i2 , nc , ncm , ncp , ncx ;
    double hx , hx2 , tau , gam , s0 , s1 , s2 , sk0 , sk1 , sk2 , sf0 ,
    sf1 , sf2 , sy1 , sy2 , sss , y0m , y0p , timestop ;
    double * xx , * aa , * bb , * cc , * ff , * y0 , * y1 , * yt , * y2 , * y3 ,
    * y4 , * al ;
    timestop = 1.0;
    MyNetInit(&argc,&argv,&np,&mp,&nl,pname,&tick);

    fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
    sleep(1);
    if ( mp ==0) 
    {
        sprintf ( sname , "%s.d" , vname ) ;
        ier = fopen_m (& Fi , sname , "rt" ) ;
        if ( ier !=0) mpierr ( " Data file not opened " ,2) ;
        fscanf ( Fi , "xa =%le\n" ,&xa ) ; // нижняя граница по x
        fscanf ( Fi , "xb =%le\n" ,&xb ) ; // верхняя граница по x
        fscanf ( Fi , "r0 =%le\n" ,&r0 ) ; // гран. условие
        fscanf ( Fi , "q0 =%le\n" ,&q0 ) ; // гран. условие
        fscanf ( Fi , "u0 =%le\n" ,&u0 ) ; // гран. условие
        fscanf ( Fi , "theta =%le\n" ,&theta ) ; // параметр в методе ньютона
        fscanf ( Fi , "tau1 =%le\n" ,&tau1 ) ; // шаг по времени
        fscanf ( Fi , "tmax =%le\n" ,&tmax ) ; // верхняя граница по времени
        fscanf ( Fi , "epst =%le\n" ,&epst ) ; //точность
        fscanf ( Fi , "nx =%d\n" ,&nx ) ; // количество шагов по x
        fscanf ( Fi , "ntp =%d\n" ,&ntp ) ; // вывод инфы каждый ntp шагов
        fscanf ( Fi , "ntm =%d\n" ,&ntm ) ; // количество узлов
        fscanf ( Fi , "lp =%d\n" ,&lp ) ;
        fclose_m (& Fi ) ;
        if ( argc >1) sscanf ( argv [1] , " % d " ,& nx ) ;
        if ( argc >2) sscanf ( argv [2] , " % d " ,& ntp ) ;
        if ( argc >3) sscanf ( argv [3] , " % d " ,& ntm ) ;
        if ( argc >4) sscanf ( argv [4] , " % lf " ,& timestop ) ;
    }
    // xa=0;
    // xb=1;
    // theta = 0.5;
    // r0=0;
    // q0=0;
    // u0=0;
    // theta=0.5;
    // tau1=0.01;
    // tmax=2.0;
    // epst=1e-9;
    // nx=400;
    // ntp=100;
    // ntm=1000000;
    // lp=0;
    if ( np >1) 
    {
        if ( mp ==0) 
        {
            buf.ddata[0] = xa ;
            buf.ddata[1] = xb ;
            buf.ddata[2] = r0 ;
            buf.ddata[3] = q0 ;
            buf.ddata[4] = u0 ;
            buf.ddata[5] = theta ;
            buf.ddata[6] = tau1 ;
            buf.ddata[7] = tmax ;
            buf.ddata[8] = epst ;
            buf.ddata[9] = nx ;
            buf.ddata[10] = ntp ;
            buf.ddata[11] = ntm ;
            buf.ddata[12] = lp ;
        }
        MPI_Bcast ( buf.ddata ,16 , MPI_DOUBLE ,0 , MPI_COMM_WORLD ) ;
        if ( mp >0) 
            {
                xa= buf.ddata[0];
                xb= buf.ddata[1];
                r0= buf.ddata[2];
                q0= buf.ddata[3];
                u0= buf.ddata[4];
                theta = buf.ddata[5];
                tau1 = buf.ddata[6];
                tmax = buf.ddata[7];
                epst = buf.ddata[8];
                nx = buf.ddata[9];
                ntp = buf.ddata[10];
                ntm = buf.ddata[11];
                lp = buf.ddata[12];
            }
    }
    fprintf ( stderr , " Netsize : %d , process : %d , system : %s , tick 7=%12le \n " , np , mp , pname , tick ) ;
    fprintf ( stderr , " xa =%le xb =%le r0 =%le \n " ,xa , xb , r0 ) ;
    fprintf ( stderr , " q0 =%le u0 =%le theta =%le \n " ,q0 , u0 , theta ) ;
    fprintf ( stderr , " tau1 =%le tmax =%le epst =%le \n " , tau1 , tmax , epst );
    fprintf ( stderr , " nx =% d ntp =% d ntm =% d lp =% d \n " ,nx , ntp , ntm , lp ) ;
    t1 = MPI_Wtime () ;
    omg1 = 1.0 / tau1 ;
    hx = ( xb - xa ) / nx ; hx2 = hx * hx ;
    tau = 0.5 * hx / sqrt ( MAX_K ) ;
    gam = tau / hx2 ;
    s0 = dmin ( tmax / tau ,1000000000.0);
    ntm = imin ( ntm ,( int ) s0 );
    fprintf ( stderr , " omg1 =%le \n " , omg1 ) ;
    fprintf ( stderr , " hx =%le tau =%le ntm =% d \n " ,hx , tau , ntm ) ;
    if ( mp == 0) fprintf ( stderr , " nx =% d hx =%le tau =%le ntm =% d \n " ,nx , hx , tau , ntm ) ;
    if ( mp ==0) mp_l = -1; else mp_l = mp - 1;
    if ( mp == np -1) mp_r = -1; else mp_r = mp + 1;
    MyRange ( np , mp ,0 , nx ,& i1 ,& i2 ,& nc ) ;
    ncm = nc -1; ncp = 2*( np -1) ; ncx = imax ( nc , ncp ) ;
    fprintf ( stderr , " i1 =% d i2 =% d nc =% d \n " ,i1 , i2 , nc ) ;
    xx = ( double *) ( malloc ( sizeof ( double ) * nc ) ) ;
    y0 = ( double *) ( malloc ( sizeof ( double ) * nc ) ) ;
    y1 = ( double *) ( malloc ( sizeof ( double ) * nc ) ) ;
    yt = ( double *) ( malloc ( sizeof ( double ) * nc ) ) ;
    aa = ( double *) ( malloc ( sizeof ( double ) * nc ) ) ;
    bb = ( double *) ( malloc ( sizeof ( double ) * nc ) ) ;
    cc = ( double *) ( malloc ( sizeof ( double ) * nc ) ) ;
    ff = ( double *) ( malloc ( sizeof ( double ) * nc ) ) ;
    al = ( double *) ( malloc ( sizeof ( double ) * ncx ) ) ;
    if ( np >1) {
        y2 = ( double *) ( malloc ( sizeof ( double ) * nc ) ) ;
        y3 = ( double *) ( malloc ( sizeof ( double ) * nc ) ) ;
        y4 = ( double *) ( malloc ( sizeof ( double ) *9* ncp ) ) ;
    }
    for ( i =0; i < nc ; i ++) xx[i] = xa + hx * ( i1 + i ) ; // grid
    for ( i =0; i < nc ; i ++) y1[i] = g0 ( xx[i]) ; // initial profile
    for ( i =0; i < nc ; i ++) yt[i] = y1[i];
    ntv = 0; tv = 0.0; gt = 1.0;
    // Time loop :
    do {
        ntv ++; tv += tau ;
        nyv = 0;
        do {
            nyv ++;
            for ( i =0; i < nc ; i ++) y0[i] = y1[i];
            if ( np >1)
            {
                BndAExch1D ( mp_l ,1 , y0 + 0 ,& y0m , mp_r ,1 , y0 + ncm ,& y0p ) ;
            }
            else
            { 
                y0m = 0.0; 
                y0p = 0.0; 
            }
            for ( i =0; i < nc ; i ++) 
            {
                ii = i1 + i ;
                if ( ii ==0) 
                {
                    aa[i] = 0.0; bb[i] = 0.0; cc[i] = 1.0;
                    ff[i] = g1 ( tv )+ tau*f(xx[ i ],tv);
                }
                else if ( ii == nx ) 
                {
                    aa[i] = 0.0; bb[i] = 0.0; cc[i] = 1.0;
                    ff[i] = g2 ( tv )+ tau*f(xx[ i ],tv) ;
                }
                else 
                {
                    s0 = k ( xx[i] , y0[i]) ;
                    sk0 = k1 ( y0[i]) ;
                    if ( i ==0) 
                    {
                        s1 = k ( xx[i] - hx , y0m ) ;
                        sk1 = k1 ( y0m ) ;
                        sy1 = y0[i] - y0m ;
                    }
                    else {
                        s1 = k ( xx[i-1] , y0[i-1]) ;
                        sk1 = k1 ( y0[i-1]) ;
                        sy1 = y0[i] - y0[i-1];
                        }
                    if ( i == ncm ) 
                    {
                        s2 = k ( xx[i]+ hx , y0p ) ;
                        sk2 = k1 ( y0p ) ;
                        sy2 = y0[i] - y0p ;
                    }
                    else 
                    {
                        s2 = k ( xx [ i +1] , y0 [ i +1]) ;
                        sk2 = k1 ( y0 [ i +1]) ;
                        sy2 = y0[i] - y0 [ i +1];
                    }
                    aa[i] = gam * 2.0 * s0 * s1 / ( s0 + s1 ) ;
                    bb[i] = gam * 2.0 * s0 * s2 / ( s0 + s2 ) ;
                    cc[i] = 1.0 + aa[i] + bb[i] + tau * q ( y0[i]) ;
                    ff[i] = y0[i] + tau*f(xx[i],tv);
                    // Добавка к  A
                    sss = s0 / ( s0 + s1 ) ;
                    sf1 = sss * sss ;
                    sf1 = 2 * gam * sk1 * sf1 * sy1 ;
                    aa[i] = aa[i] - theta * sf1 ;
                    // Добавка к B
                    sss = s0 / ( s0 + s2 ) ;
                    sf2 = sss * sss ;
                    sf2 = 2 * gam * sk2 * sf2 * sy2 ;
                    bb[i] = bb[i] - theta * sf2 ;
                    // Добавка к  C
                    sss = s1 / ( s0 + s1 ) ;
                    sf0 = sss * sss * sy1 ;
                    sss = s2 / ( s0 + s2 ) ;
                    sf0 += sss * sss * sy2 ;
                    sf0 *= 2 * gam * sk0 ;
                    sf0 += tau * q1 ( y0[i]) - 1;
                    cc[i] = cc[i] + theta * sf0 ;
                    // addition to F
                    if ( i ==0) 
                    {
                        sss = sf1 * y0m + sf0* y0[i] + sf2 * y0 [ i +1];
                    }
                    else if ( i == ncm ) 
                    {
                        sss = sf1 * y0[i-1] + sf0 * y0[i] + sf2 * y0p ;
                    }
                    else 
                    {
                        sss = sf1 * y0[i-1] + sf0 * y0[i] + sf2 * y0 [ i +1];
                    }
                    ff[i] = ff[i] + theta * sss ;
                }
            }
            if ( np <2) 
            {
                ier = prog_right ( nc , aa , bb , cc , ff , al , y1 ) ;
            }
            else
            {
                ier = prog_rightpm ( np , mp , nc , ntv , aa , bb , cc , ff, al , y1 , y2 , y3 , y4 ) ;
            }
            if ( ier !=0) mpierr ( " Bad solution " ,1) ;
            gt = 0.0;
            for ( i =0; i < nc ; i ++) 
            {
                if ( y0[i]!=0) { s0 = ( y1[i]/ y0[i] -1.0) ; }
                else if ( y1[i]==0) { s0 = 0.0; }
                else { s0 = 1.0; }
                // s0 = ( y1[i] - y0[i]) ;
                gt = dmax ( gt , dabs ( s0 ) ) ;
            }
            gt = gt / tau ;
            if ( np >1) 
            {
                s0 = gt ; MPI_Allreduce (& s0 ,& gt ,1 , MPI_DOUBLE , MPI_MAX, MPI_COMM_WORLD ) ;
            }
        } while (( nyv <1000) && ( gt > epst ) ) ;
        if ( ntv % ntp == 0) 
        {
            if ( mp == 0) 
            {
                t2 = MPI_Wtime () - t1 ;
                fprintf ( stderr , " ntv =%d tv =%le gt =%le nyv =%4d tcpu =%le \n " ,ntv , tv , gt , nyv , t2 ) ;
            }
        }
        for ( i =0; i < nc ; i ++){ yt[i] = y1[i];}
        if ( lp >0) 
        {
            fprintf ( stderr , " ntv =%d tv =%le gt =%le \n " ,ntv , tv , gt ) ;
            for ( i =0; i < nc ; i ++)
                fprintf ( stderr , " i =%8d x =%12le y1 =%12le \n " ,( i1 + i ) , xx[i] ,yt[i]) ;
        }
    } while (( ntv < ntm ) && ( tv < timestop ) ) ;
    t1 = MPI_Wtime () - t1 ;
    sprintf ( sname , " % s_ %02 d . dat " , vname , np ) ;
    OutFun1DP ( sname , np , mp , nc , xx , yt ) ;
    fprintf ( stderr , " ntv =% d tv =%le gt =%le time =%le \n " ,ntv , tv , gt , t1);
    if ( mp == 0) fprintf ( stderr , " ntv =% d tv =%le gt =%le tcpu =%lf \n " ,ntv , tv , gt , t1 ) ;
    // ier = fclose_m (& Fo ) ;
    MPI_Finalize () ;
    return 0;
}
