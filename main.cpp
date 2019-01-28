#include "def.h"
#include <stdlib.h>
#include <iostream>
#include <string>
//#include <time.h>
using namespace std;
#ifndef max
#define max(a,b) (( (a) >= (b)) ? (a) : (b))
#endif
#define N 16 //8个变量

void equations(double var[],double diff[],double spin,double krz_d[],double massratio);
void dboydcar(double spin, double car[], double boy[], double dboy_dcar[][4]);
void dboydcar_boyknown(double spin, double car[], double boy[], double dboy_dcar[][4]);
void dcardboy(double spin, double boy[], double car[], double dcar_dboy[][4]);
void highderiv(double spin, double krz_d[], double boy[], double boy_taudot[], double x[][4]);
void radacc(double spin, double x[][4], double acc[]);

int main(int argc, char *argv[])
{//RK45的参数
	double a1 = 1.0 / 4.0;
	double b1 = 3.0 / 32.0;
	double b2 = 9.0 / 32.0;
	double c1 = 1932.0 / 2197.0;
	double c2 = -7200.0 / 2197.0;
	double c3 = 7296.0 / 2197.0;
	double d1 = 439.0 / 216.0;
	double d2 = -8.0;
	double d3 = 3680.0 / 513.0;
	double d4 = -845.0 / 4104.0;
	double e1 = -8.0 / 27.0;
	double e2 = 2.0;
	double e3 = -3544.0 / 2565.0;
	double e4 = 1859.0 / 4104.0;
	double e5 = -11.0 / 40.0;
	double x1 = 25.0 / 216.0;
	double x2 = 0.0;
	double x3 = 1408.0 / 2565.0;
	double x4 = 2197.0 / 4104.0;
	double x5 = -1.0 / 5.0;
	double z1 = 16.0 / 135.0;
	double z2 = 0.0;
	double z3 = 6656.0 / 12825.0;
	double z4 = 28561.0 / 56430.0;
	double z5 = -9.0 / 50.0;
	double z6 = 2.0 / 55.0;

	double spin, spin2, krz_d[10] = { 0 };
	double isco, xin;
	double robs_i, robs_f;
	double pstep;
	int i, j, k, m;
	int ii,jj;
	char filename_o[128];

	FILE *foutput, *finput;

	double r, th, t, phi, ur, ut, uth, uphi, tau=0, dtau; //坐标，4速度，固有时，固有时步长

	//double rnew, thnew, tnew, phinew, urnew, utnew, uthnew, uphinew, dtau;

	double E, Lz;

	double F_theta, F_r,F_t,F_phi;

	double g[4][4] = { 0 };
	double Gamma[4][4][4] = { 0 };
	double u[4] = { 0 };
	double k1[N], k2[N], k3[N], k4[N],k5[N],k6[N];//16个变量的顺序是t,r,\theta,\phi,u^t, u^r,u^\theta, u^\phi, ts,rs,\theta s,\phi s,us^t, us^r,us^\theta, us^\phi
	double var[N];
	double ita,mu=-1;//计算过程中四速度的模是\ita， 描述粒子性质的是mu 有质量：mu=-1，无质量： mu=0

//input argument: 1. M, 2. spin, 3. E, 4. Lz, 5. Q, 6. r0, 7. tottime, 8-10. krz_d, 11. mass ratio
	
	/***************************改成秒间隔需要改的部分********************/

	double dt;
	double Grav, clight, Msol, M, dt_sec = 0.1;
	double t_sec = 0;
	
	Grav = 6.674e-11;//#引力常数
	clight = 2.998e8;//#光速
	Msol = 1.989e30;  //#太阳质量，以千克做单位
	
	M = atof(argv[1]);
	//dt_sec = dt*M*Msol*Grav / clight / clight / clight

	dt = dt_sec*clight*clight*clight / M / Msol / Grav;
	
	/***************************改成秒间隔需要改的部分********************/

	dtau = dt;//初始步长
	double tottime=atof(argv[7]);
	double massratio = atof(argv[11]);
	spin = atof(argv[2]);//+0.001*0.5;
	//for (spin = 0.55428452790227312;spin <  0.55428452790227312 +0.01;spin += 0.05) {
		printf("spin=%f \n", spin);
		krz_d[1] = atof(argv[8]),krz_d[2]=atof(argv[9]), krz_d[3]=atof(argv[10]);//+0.001*0.2;
//		for (krz_d[2] = 0.0; krz_d[2] <= 0.0; krz_d[2] += 0.1) {
			printf("delta: %f, %f, %f\n",krz_d[1], krz_d[2],krz_d[3]);
			t = 0;
			//r = 13.0;
			th = Piby2;
			phi = 0;


			//ut = ( E*g[3][3]+Lz*g[0][3] ) / ( g[0][3]*g[0][3] - g[0][0]*g[3][3] );
			//uphi = (E*g[0][3] + Lz*g[0][0]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]) ;
			//uth = sqrt((-1 - g[0][0] * ut*ut - 2 * g[0][3] * ut*uphi - g[3][3] * uphi*uphi - g[1][1] * ur*ur) / (g[2][2]));

			/********************↓圆轨道的能量和角动量 from https://arxiv.org/pdf/1105.2959.pdf （好像是错的。。）↓******************/
			/*	Lz = fabs(spin*spin + 2 * spin*sqrt(r) + r*r) / sqrt(r*r*(r - 3) - 2 * spin*sqrt(r*r*r));//Kerr_circular_orbit, corotating
			E = ((pow(r, 1.25)*fabs((spin*spin + r*(r - 2))*(spin - sqrt(r*r*r))) / sqrt((r - 3)*sqrt(r) - 2 * spin)) + 2 * spin*r*Lz) / (r*(r*r*r + spin*spin*(r + 2))); //Kerr_circular_orbit, corotating

			Lz = fabs(spin*spin - 2 * spin*sqrt(r) + r*r) / sqrt(r*r*(r - 3) + 2 * spin*sqrt(r*r*r));//Kerr_circular_orbit, counterrotating
			E = ((pow(r, 1.25)*fabs((spin*spin + r*(r - 2))*(spin + sqrt(r*r*r))) / sqrt((r - 3)*sqrt(r) + 2 * spin)) + 2 * spin*r*Lz) / (r*(r*r*r + spin*spin*(r + 2))); //Kerr_circular_orbit, counterrotating
			*/
			/********************↑from https://arxiv.org/pdf/1105.2959.pdf （好像是错的。。。）↑******************/

			/********************↓圆轨道的能量和角动量 from http://www.tapir.caltech.edu/~chirata/ph236/2011-12/lec27.pdf （这个是对的）↓******************/
			//Lz = ( sqrt(r) -2*spin/r + spin*spin/sqrt(r*r*r) ) / sqrt(1 - 3 / r + 2 * spin / sqrt(r*r*r));
			//E = (1-2/r + spin /sqrt(r*r*r) ) / sqrt(1-3/r + 2*spin/sqrt(r*r*r) );
			/********************↑from http://www.tapir.caltech.edu/~chirata/ph236/2011-12/lec27.pdf ↑******************/


			double horizon, r0, th0 = Piby2, phi0 = 0, t0 = 0;//初始条件
			double ecc, p;//eccentricity和rectum 
			double rmax, rmin, invgmax[4][4], invgmin[4][4];//r的上下限以及该处的g
			double invg[4][4];
			double EoverL;//由e和p算E和Lz的中间变量
			double EoverL2, E2, L2;
			//double iota = Pi / 6;
			double Q,Q0;//carter constant
			horizon = 1 + sqrt(1 - spin*spin);
			/***************************赤道面上由e,p决定的轨道↓******************/
			find_isco_KRZ(spin, krz_d, 100., isco);

			//for (ecc = 0.41149551640029858;ecc <= 0.41149551640029858;ecc = ecc + 0.004) {
			/*for (ecc = 0.5;ecc <= 0.5;ecc = ecc + 0.004) {
				//	for (p = 6.4825501282607396;p <= 6.4825501282607396;p = p + 0.04) {
				for (p = 6;p <= 6;p = p + 0.04) {*/
					
			//ecc = 0.43173473149300562,p= 7.1717260434110983;
					E = atof(argv[3]), Lz=atof(argv[4]), Q=atof(argv[5]);//python解出来的   and note that it's not useful actually, Q is calculated below, and for non-Kerr case this is just initial Q
					th0 = Piby2;//赤道面上carter constant和theta方向速度关系比较简单
					r0 = atof(argv[6]);//看看取初始在rmax会怎么样？
					Q0 = Q;
			/*
					rmax = p / (1 - ecc);
					rmin = p / (1 + ecc);
					r0 = rmax;
					metric_KRZ_inverse(spin, krz_d, rmax, th0, invgmax);
					metric_KRZ_inverse(spin, krz_d, rmin, th0, invgmin);


					EoverL = ((invgmax[3][0] - invgmin[3][0]) + sqrt((invgmax[3][0] - invgmin[3][0]) *(invgmax[3][0] - invgmin[3][0]) - (invgmax[0][0] - invgmin[0][0])*(invgmax[3][3] - invgmin[3][3]))) / (invgmax[0][0] - invgmin[0][0]);
					Lz = sqrt((invgmax[3][0] - invgmin[3][0]) / (EoverL*EoverL*(invgmin[3][0] * invgmax[0][0] - invgmax[3][0] * invgmin[0][0]) + (invgmin[3][0] * invgmax[3][3] - invgmax[3][0] * invgmin[3][3])));

					E = EoverL*Lz;

					*/
					//EoverL2 = ((invgmax[3][0] - invgmin[3][0]) - sqrt((invgmax[3][0] - invgmin[3][0]) *(invgmax[3][0] - invgmin[3][0]) - (invgmax[0][0] - invgmin[0][0])*(invgmax[3][3] - invgmin[3][3]))) / (invgmax[0][0] - invgmin[0][0]);
					//L2 = sqrt((invgmax[3][0] - invgmin[3][0]) / (EoverL*EoverL*(invgmin[3][0] * invgmax[0][0] - invgmax[3][0] * invgmin[0][0]) + (invgmin[3][0] * invgmax[3][3] - invgmax[3][0] * invgmin[3][3])));
					//
					//E2 = abs(EoverL2*Lz);
					//if (E2 < E) {
					//	E = E2;
					//}
					r = r0;
					th = th0;
					t = t0;
					phi = phi0;
					metric_KRZ(spin, krz_d, r, th, g);
					metric_KRZ_inverse(spin, krz_d, r, th, invg);



					printf("template E=%.6f Lz=%.6f Q=%.6f\n",  E, Lz, Q);
					/********************↓二分找E,Lz↓******************/
					/*
					if (abs(ecc - 0.0) < 1e-6) {
						r = p;r0 = p;
						double upE = 1, downE = 0.9, eps = 1e-10;
						double invg[4][4];
						double curf, upf, downf;
						Christoffel_KRZ(spin, krz_d, r, th, Gamma);
						metric_KRZ_inverse(spin, krz_d, r, th, invg);

						E = upE;
						Lz = (invg[0][3] * E + sqrt((invg[0][3] * E)*(invg[0][3] * E) - invg[3][3] * (invg[0][0] * E*E + 1))) / invg[3][3];

						ut = (E*g[3][3] + Lz*g[0][3]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
						uphi = -(E*g[0][3] + Lz*g[0][0]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
						upf = Gamma[1][0][0] * ut*ut + 2 * Gamma[1][0][3] * ut*uphi + Gamma[1][3][3] * uphi*uphi;


						E = downE;
						Lz = (invg[0][3] * E + sqrt((invg[0][3] * E)*(invg[0][3] * E) - invg[3][3] * (invg[0][0] * E*E + 1))) / invg[3][3];

						ut = (E*g[3][3] + Lz*g[0][3]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
						uphi = -(E*g[0][3] + Lz*g[0][0]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
						downf = Gamma[1][0][0] * ut*ut + 2 * Gamma[1][0][3] * ut*uphi + Gamma[1][3][3] * uphi*uphi;


						for (;;) {
							E = 0.5*(upE + downE);
							Lz = (invg[0][3] * E + sqrt((invg[0][3] * E)*(invg[0][3] * E) - invg[3][3] * (invg[0][0] * E*E + 1))) / invg[3][3];

							ut = (E*g[3][3] + Lz*g[0][3]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
							uphi = -(E*g[0][3] + Lz*g[0][0]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
							curf = Gamma[1][0][0] * ut*ut + 2 * Gamma[1][0][3] * ut*uphi + Gamma[1][3][3] * uphi*uphi;

							if (curf*downf < 0) {
								upE = E;
								upf = curf;
							}
							else {
								downE = E;
								downf = curf;
							}
							if (abs(upE - downE) < eps) {
								E = 0.5*(upE + downE);
								Lz = (invg[0][3] * E + sqrt((invg[0][3] * E)*(invg[0][3] * E) - invg[3][3] * (invg[0][0] * E*E + 1))) / invg[3][3];
								break;
							}

						}

					}*/

					/********************↑二分找圆轨道 ↑******************/
					/*
					double testE, testL, testut, testup, utoverup;
					if (abs(ecc - 0.0) < 1e-6) {
						///解方程找圆轨道↓
						//p = 10;
						Christoffel_KRZ(spin, krz_d, p, th, Gamma);
						metric_KRZ(spin, krz_d, p, th, g);
						utoverup = (-Gamma[1][0][3] + sqrt(Gamma[1][0][3] * Gamma[1][0][3] - Gamma[1][3][3] * Gamma[1][0][0])) / Gamma[1][0][0];
						testup = sqrt(-1 / (utoverup*utoverup*g[0][0] + 2 * utoverup*g[0][3] + g[3][3]));
						testut = testup*utoverup;
						testE = -g[0][0] * testut - g[0][3] * testup;
						testL = g[0][3] * testut + g[3][3] * testup;
						E = testE;
						Lz = testL;
						//   解方程找圆轨道↑
					}
					*/

					/*looking for circular orbit↓*/
					/*
					metric_KRZ(spin, krz_d, r, th, g);
					Christoffel_KRZ(spin, krz_d, r, th, Gamma);
					sprintf(filename_o, "test.dat");
					foutput = fopen(filename_o, "w");

					for (E = 0.5;E < 1;E += 0.001) {
						for (Lz = 0;Lz < 10;Lz += 0.01) {
							ut = (E*g[3][3] + Lz*g[0][3]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
							uphi = -(E*g[0][3] + Lz*g[0][0]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
							fprintf(foutput, "%.6f \t ", Gamma[1][0][0] * ut*ut + 2 * Gamma[1][0][3] * ut*uphi + Gamma[1][3][3] * uphi*uphi);

						}
						for (Lz = 0;Lz < 10;Lz += 0.01) {
							ut = (E*g[3][3] + Lz*g[0][3]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
							uphi = -(E*g[0][3] + Lz*g[0][0]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
							fprintf(foutput, "%.6f \t ", g[0][0] * ut*ut + 2 *g[0][3] * ut*uphi + g[3][3] * uphi*uphi);

						}
						fprintf(foutput, "\n");
						printf("%.2f %.2f \n", E, Lz);

					}
					fclose(foutput);
					abort();

					/*looking for circular orbit↑*/

					ut = (E*g[3][3] + Lz*g[0][3]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
					uphi = -(E*g[0][3] + Lz*g[0][0]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
					//uth = sqrt(Q) / g[2][2];
					//ur = sqrt((mu - g[0][0] * ut*ut - 2 * g[0][3] * ut*uphi - g[3][3] * uphi*uphi - g[2][2] * uth*uth) / (g[1][1]));
					ur = 0;
					if (abs((mu - g[0][0] * ut*ut - 2 * g[0][3] * ut*uphi - g[3][3] * uphi*uphi - g[1][1] * ur*ur) / (g[2][2])) < eps) uth = 0;
					else uth = sqrt((mu - g[0][0] * ut*ut - 2 * g[0][3] * ut*uphi - g[3][3] * uphi*uphi - g[1][1] * ur*ur) / (g[2][2]));
					//uth = 0;
					//ur = 0;
					Q = uth*uth*g[2][2] * g[2][2];
					E = -g[0][0] * ut - g[0][3] * uphi;
					Lz = g[0][3] * ut + g[3][3] * uphi;
					u[0] = ut;
					u[1] = ur;
					u[2] = uth;
					u[3] = uphi;

					//printf("template E=%.6f Lz=%.6f\n", ecc, p, E, Lz);

					//sprintf(filename_o, "circular_trace_spin%.2f_d%.2f_r%.2f.dat", spin, krz_d,r);

					printf("template E=%.6f Lz=%.6f Q=%.6f\n", E, Lz, Q);

					//测试坐标变换
					/*
					double car[4]; double boy[4];
					double dboy_dcar[4][4]; double dcar_dboy[4][4];

					boy[0] = t;boy[1] = r;boy[2] = th;boy[3] = phi;
					dcardboy(spin, boy, car, dcar_dboy);
					cout << car[0]<<car[1]<<car[2]<<car[3];

					dboydcar(spin, car, boy, dboy_dcar);
					int flag = 0;
					for (int ind = 0;ind < 4;ind++) {
						if (!isfinite(boy[ind])) {
							flag = 1;
							cout << "boy [" << ind << ']' << endl;
						}
						if (!isfinite(car[ind])) {
							flag = 1;
							cout << "car [" << ind << ']' << endl;
						}
						for (int ind2 = 0;ind2 < 4;ind2++) {
							if (!isfinite(dboy_dcar[ind][ind2])) {
								flag = 1;
								cout << "dboy_dcar [" << ind << "]["<<ind2<<']' << endl;
							}
							if (!isfinite(dcar_dboy[ind][ind2])) {
								flag = 1;
								cout << "dcar_dboy [" << ind << "]["<<ind2<<']' << endl;
							}
						}
					}
					if (flag == 1) system("pause");*/
					//上面在测试坐标变换

					sprintf(filename_o, "trace_M%.0f_spin%.6f_E%.6f_Lz%.6f_Q%.6f_d1%.6f_d2%.6f_d3%.6f.dat", M, spin, E, Lz, Q, krz_d[1], krz_d[2], krz_d[3]);


					foutput = fopen(filename_o, "w");
					//Zprintf("index\t tau\t t\t r\t theta\t phi\t ut\t ur\t uth\t uphi\t ita\t E\t Lz\n");

					var[0] = t;	var[1] = r;	var[2] = th;	var[3] = phi;	var[4] = ut;	var[5] = ur;	var[6] = uth;	var[7] = uphi;
					var[8] = 0;	var[9] = 0;	var[10] = 0;	var[11] = 0;	var[12] = 0;	var[13] = 0;	var[14] = 0;	var[15] = 0;//小量初始为0
					double E0 = E, L0 = Lz;
					double h;//RK45, adaptive step
					double diff[N];
					double vars_temp[N];//temp used in RK45
					double z[N], y[N];//used to estimate error
					double err, maxerr=1e-10, minerr=1e-11;//计算过程中的误差，容许的最大误差，容许的最小误差
					int check=0;//标记有没有超过最大误差或小于最小误差
					int index = 0;
					h = dtau;
					tau = 0;
					//time_t optime = time(NULL);
					for(;t_sec<tottime;){


						t_sec = var[0] * M*Msol*Grav / clight / clight / clight;

						metric_KRZ(spin, krz_d, var[1]+var[9], var[2]+var[10], g);
						ut = var[4]+var[12];ur = var[5]+var[13];uth = var[6]+var[14];uphi = var[7]+var[15];
						
						ita = g[0][0] * ut*ut + g[1][1] * ur*ur + g[2][2] * uth*uth + g[3][3] * uphi*uphi + 2 * g[3][0] * uphi*ut;
						Q = uth*g[2][2] * uth*g[2][2] + cos(var[2])*cos(var[2])*(spin*spin*(ita*ita - E*E) + Lz*Lz / sin(var[2]) / sin(var[2]));
						E = -g[0][0] * ut - g[0][3] * uphi;
						Lz = g[0][3] * ut + g[3][3] * uphi;
						
						//下面在测试坐标变换
						/*
						double car[4]; double boy[4];
						double dboy_dcar[4][4]; double dcar_dboy[4][4];
						boy[0] = var[0];boy[1] = var[1];boy[2] = var[2];boy[3] = car[3];
						dcardboy(spin, boy, car, dcar_dboy);
						cout << "car:   "<< car[0] <<'\t'<< car[1] <<'\t'<< car[2] <<'\t'<< car[3]<<endl;

						dboydcar(spin, car, boy, dboy_dcar);
						int flag = 0;
						for (int ind = 0;ind < 4;ind++) {
							if (!isfinite(boy[ind])) {
								flag = 1;
								cout << "boy [" << ind << ']' << endl;
							}
							if (!isfinite(car[ind])) {
								flag = 1;
								cout << "car [" << ind << ']' << endl;
							}
							for (int ind2 = 0;ind2 < 4;ind2++) {
								if (!isfinite(dboy_dcar[ind][ind2])) {
									flag = 1;
									cout << "dboy_dcar [" << ind << "][" << ind2 << ']' << endl;
								}
								if (!isfinite(dcar_dboy[ind][ind2])) {
									flag = 1;
									cout << "dcar_dboy [" << ind << "][" << ind2 << ']' << endl;
								}
							}
						}
						if (flag == 1) system("pause");*/
						//上面在测试坐标变换

						check = 0;
						equations(var, diff,spin,krz_d,massratio);
						for (i = 0; i < N; i++)
						{
							k1[i] = h*diff[i];
							vars_temp[i] = var[i] + a1*k1[i];
						}

						equations(vars_temp, diff,spin,krz_d,massratio);
						for (i = 0; i < N; i++)
						{
							k2[i] = h*diff[i];
							vars_temp[i] = var[i] + b1*k1[i] + b2*k2[i];
						}

						equations(vars_temp, diff,spin,krz_d, massratio);
						for (i = 0; i < N; i++)
						{
							k3[i] = h*diff[i];
							vars_temp[i] = var[i] + c1*k1[i] + c2*k2[i] + c3*k3[i];
						}

						equations(vars_temp, diff,spin,krz_d, massratio);
						for (i = 0; i < N; i++)
						{
							k4[i] = h*diff[i];
							vars_temp[i] = var[i] + d1*k1[i] + d2*k2[i] + d3*k3[i] + d4*k4[i];
						}

						equations(vars_temp, diff,spin,krz_d, massratio);
						for (i = 0; i < N; i++)
						{
							k5[i] = h*diff[i];
							vars_temp[i] = var[i] + e1*k1[i] + e2*k2[i] + e3*k3[i] + e4*k4[i] + e5*k5[i];
						}

						equations(vars_temp, diff,spin,krz_d, massratio);
						for (i = 0; i < N; i++)
							k6[i] = h*diff[i];

						for (i = 0; i < N; i++)
						{
							y[i] = var[i] + x1*k1[i] + x2*k2[i] + x3*k3[i] + x4*k4[i] + x5*k5[i];
							z[i] = var[i] + z1*k1[i] + z2*k2[i] + z3*k3[i] + z4*k4[i] + z5*k5[i] + z6*k6[i];
							err = fabs((y[i] - z[i]) / max(fabs(var[i]), fabs(y[i]) ) );
							if (err > maxerr) {
								check = 1;//有些误差太大
							}
							else if (err < minerr&&check != 1) {
								check = -1;//所有误差都太小
							}
						}

						if (check == 1) {
							h /= 1.1;
						}
						else if (check == -1) {
							if (check != 1) {
								fprintf(foutput, "%d\t %.10f\t %.10f\t %.10f\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e \n", index, t_sec, tau, var[0], var[1], var[2], var[3], var[4], var[5], var[6], var[7], /*F_t*/(z1*k1[4] + z2*k2[4] + z3*k3[4] + z4*k4[4] + z5*k5[4] + z6*k6[4]), /*F_r*/(z1*k1[5] + z2*k2[5] + z3*k3[5] + z4*k4[5] + z5*k5[5] + z6*k6[5]), /*F_theta*/z1*k1[6] + z2*k2[6] + z3*k3[6] + z4*k4[6] + z5*k5[6] + z6*k6[6], /*F_phi*/z1*k1[7] + z2*k2[7] + z3*k3[7] + z4*k4[7] + z5*k5[7] + z6*k6[7]);
								printf("%d\t %.10f\t %.3f\t %.3f\t %.6f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.10e \t %.10e \t %.10e \t %.10e\n", index, t_sec, tau, var[0], var[1], var[2], var[3], var[4]/*ut*/, var[5]/*ur*/, var[6]/*uth*/, var[7]/*uphi*/, ita - mu, (E - E0) / E0, (Lz - L0) / L0, Q);
								fflush(foutput);
							}
							tau = tau + h;
							h *= 1.1;
							index++;
							for (i = 0; i < N; i++)
								var[i] = y[i];
						}
						else {
							if (check != 1) {
								fprintf(foutput, "%d\t %.10f\t %.10f\t %.10f\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e \n", index, t_sec, tau, var[0], var[1], var[2], var[3], var[4], var[5], var[6], var[7], /*F_t*/(z1*k1[4] + z2*k2[4] + z3*k3[4] + z4*k4[4] + z5*k5[4] + z6*k6[4]), /*F_r*/(z1*k1[5] + z2*k2[5] + z3*k3[5] + z4*k4[5] + z5*k5[5] + z6*k6[5]), /*F_theta*/z1*k1[6] + z2*k2[6] + z3*k3[6] + z4*k4[6] + z5*k5[6] + z6*k6[6], /*F_phi*/z1*k1[7] + z2*k2[7] + z3*k3[7] + z4*k4[7] + z5*k5[7] + z6*k6[7]);
								//printf("%d\t %.10f\t %.3f\t %.3f\t %.6f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.10e \t %.10e \t %.10e \t %.10e\n", index, t_sec, tau, var[0], var[1], var[2], var[3], var[4]/*ut*/, var[5]/*ur*/, var[6]/*uth*/, var[7]/*uphi*/, ita - mu, (E - E0) / E0, (Lz - L0) / L0, Q);
								fflush(foutput);
							}
							tau = tau + h;
							index++;
							for (i = 0; i < N; i++)
								var[i] = y[i];
						}
						

					}
					fclose(foutput);
					//time_t edtime = time(NULL);
					//cout << (edtime - optime) << endl;
					//system("pause");
				
			
		//}

	//}
	return 0;

}

void equations(double var[], double diff[],double spin,double krz_d[],double massratio) {//由8个变量算斜率
	diff[0] = var[4];diff[1] = var[5];diff[2] = var[6];diff[3] = var[7];
	double F_r, F_t, F_theta, F_phi;
	int ii, jj;
	double Gamma[4][4][4];
	double u[4]; 
	u[0] = var[4]+var[12];u[1] = var[5]+var[13];u[2] = var[6]+var[14];u[3] = var[7]+var[15];

	Christoffel_KRZ(spin, krz_d, var[1], var[2], Gamma);
	F_r = 0;
	F_theta = 0;
	F_t = 0;
	F_phi = 0;
	for (ii = 0;ii < 4;ii++) {
		for (jj = 0;jj < 4;jj++) {

			F_t -= 0.5*Gamma[0][ii][jj] * u[ii] * u[jj];
			F_r -= 0.5*Gamma[1][ii][jj] * u[ii] * u[jj];
			F_theta -= 0.5*Gamma[2][ii][jj] * u[ii] * u[jj];
			F_phi -= 0.5*Gamma[3][ii][jj] * u[ii] * u[jj];
		}
	}
	diff[4] = F_t;
	diff[5] = F_r;
	diff[6] = F_theta;
	diff[7] = F_phi;
	//检验highderiv
	/*
	double x[8][4] = { 0 }, boy[4];
	boy[0] = var[0]+var[8];boy[1] = var[1]+var[9]; boy[2] = var[2]+var[10]; boy[3] = var[3]+var[11];
	highderiv(spin, krz_d, boy, u, x);
	int flag = 0;
	for (int ind = 0;ind < 8;ind++) {
		
		for (int ind2 = 0;ind2 < 4;ind2++) {
			if (!isfinite(x[ind][ind2])) {
				flag = 1;
				cout << "x [" << ind << "][" << ind2 << ']' << endl;
			}
			
		}
	}
	if (flag == 1) system("pause");
	*/
	//上面在检验highderiv

	diff[8] = var[12];
	diff[9] = var[13];
	diff[10] = var[14];
	diff[11] = var[15];

	double x[8][4] = { 0 }, boy[4], acc[4] = { 0 }, acc_car[4] = { 0 };
	boy[0] = var[0] + var[8];boy[1] = var[1] + var[9]; boy[2] = var[2] + var[10]; boy[3] = var[3] + var[11];
	highderiv(spin, krz_d, boy, u, x);
	radacc(spin, x, acc_car);
	for (int i = 1;i < 4;i++) acc_car[i] = acc_car[i] * massratio*u[0] * u[0] - var[i + 12] * F_t/u[0];//乘上质量比才是真*加速度，另外要转换成对tau的加速度
	double car[4] = { 0 }, dboy_dcar[4][4] = { 0 };
	car[0] = boy[0];
	car[1] = sqrt(boy[1] * boy[1] + spin*spin) * sin(boy[2]) * cos(boy[3]);
	car[2] = sqrt(boy[1] * boy[1] + spin*spin) * sin(boy[2]) * sin(boy[3]);
	car[3] = boy[1] * cos(boy[2]);
	dboydcar_boyknown(spin, car, boy, dboy_dcar);
	for (int mu = 0;mu < 4;mu++) {
		for (int al = 0;al < 4;al++) {
			acc[mu]+= acc_car[al] * dboy_dcar[mu][al];
		}
	}
	diff[12] = acc[0];
	diff[13] = acc[1];
	diff[14] = acc[2];
	diff[15] = acc[3];
}

void dcardboy(double spin, double boy[], double car[], double dcar_dboy[][4]) {//从Boyer-Lindquist坐标换到Cartesian坐标
	car[0] = boy[0];
	car[1] = sqrt(boy[1] * boy[1] + spin*spin) * sin(boy[2]) *cos(boy[3]);
	car[2] = sqrt(boy[1] * boy[1] + spin*spin) * sin(boy[2]) *sin(boy[3]);
	car[3] = boy[1] * cos(boy[2]);
	dcar_dboy[0][0] = 1;
	dcar_dboy[0][1] = 0;
	dcar_dboy[0][2] = 0;
	dcar_dboy[0][3] = 0;
	dcar_dboy[1][0] = 0;
	dcar_dboy[2][0] = 0;
	dcar_dboy[3][0] = 0;
	dcar_dboy[1][1] = boy[1] / sqrt(boy[1] * boy[1] + spin*spin) * sin(boy[2]) * cos(boy[3]);//dx/dr
																							 //好像不太方便全用boy[]
	double r = boy[1], th = boy[2], phi = boy[3];
	double a = spin;
	dcar_dboy[1][2] = sqrt(r*r + spin*spin)*cos(th)*cos(phi);//dx/dth
	dcar_dboy[1][3] = -sqrt(r*r + spin*spin)*sin(th)*sin(phi);// dx/dphi
	dcar_dboy[2][1] = boy[1] / sqrt(boy[1] * boy[1] + spin*spin) * sin(boy[2]) * sin(boy[3]);//dy/dr
	dcar_dboy[2][2] = sqrt(r*r + a*a)*cos(th)*sin(phi);// dy/dth
	dcar_dboy[2][3] = sqrt(r*r + a*a)*sin(th)*cos(phi);// dy/dphi
	dcar_dboy[3][1] = cos(th);//dz/dr
	dcar_dboy[3][2] = -r*sin(th);//dz/dth
	dcar_dboy[3][3] = 0;
}
void dboydcar(double spin, double car[], double boy[], double dboy_dcar[][4]) {//从Cartesian坐标换到Boyer-Lindquist坐标
	double x = car[1],y = car[2],z = car[3];
	double a = spin;
	double r = sqrt(sqrt((-a *a + x *x + y *z + z *z)*(-a *a + x *x + y *z + z *z) + 4 * a *a * z*z) - a*a + x *x + y*y + z*z) / sqrt(2);
	double th;
	if (a == 0) {
		th = acos(z / sqrt(x*x + y*y + z*z));
	}
	else {
		double sqcth = (1 - x*x / a/a - y*y / a/a - z*z / a/a + sqrt(4 * a*a*z*z + (x*x + y*y + z*z - a*a)*(x*x + y*y + z*z - a*a)) / a / a) / sqrt(2);
		if (abs(sqcth) < eps * 10) th = Piby2;
		else {
			if (z > 0) {
				double cth = sqrt(abs(sqcth));
				th = acos(cth);
			}
			else {
				double cth = -sqrt(abs(sqcth));
				th = acos(cth);
			}
		}
	}
	double phi;
	if (x == 0) {
		phi = Piby2;
	}
	else {
		phi = atan(y / x);
	}

	boy[0] = car[0];boy[1] = r;boy[2] = th;boy[3] = phi;
	dboy_dcar[0][0] = 1;
	dboy_dcar[0][1] = 0;
	dboy_dcar[0][2] = 0;
	dboy_dcar[0][3] = 0;
	dboy_dcar[1][0] = 0;
	dboy_dcar[2][0] = 0;
	dboy_dcar[3][0] = 0;
	dboy_dcar[1][1] = r*r*r*x/( r*r*r*r+a*a*z*z );//dr/dx
	dboy_dcar[1][2] = r*r*r*y / (r*r*r*r + a*a*z*z);//dr/dy
	dboy_dcar[1][3] = r*(a*a + r*r)*z / (r*r*r*r + a*a*z*z);//dr/dz
	
	if (abs(th - Piby2) < eps) { //th趋于pi/2的极限
		dboy_dcar[2][1] = 0;
		dboy_dcar[2][2] = 0;
		dboy_dcar[2][3] = - 1 / r;
	}
	else if (abs(th - 0) < eps) { //th趋于0的极限
		dboy_dcar[2][1] = cos(phi) / sqrt(r*r + a*a);
		dboy_dcar[2][2] = sin(phi) / sqrt(r*r + a*a);
		dboy_dcar[2][3] = 0;
	}
	else {
		dboy_dcar[2][1] = x*cos(th)*cos(th) / (a*a*cos(th)*cos(th)*cos(th)*sin(th) + z*z*tan(th));
		dboy_dcar[2][2] = y*cos(th)*cos(th) / (a*a*cos(th)*cos(th)*cos(th)*sin(th) + z*z*tan(th));
		dboy_dcar[2][3] = -z *sin(th)*sin(th) / (a*a*cos(th)*cos(th)*cos(th)*sin(th) + z*z*tan(th));

	}

	if (abs(x - 0) < eps) {//x趋于0的极限
		dboy_dcar[3][1] = -1 / y;
		dboy_dcar[3][2] = 0;
	}
	else if (abs(y-0)<eps){//y趋于0的极限
		dboy_dcar[3][1] = 0;
		dboy_dcar[3][2] = 1 / x;
	}
	else {
		dboy_dcar[3][1] = -y / (x*x + y*y);
		dboy_dcar[3][2] = x / (x*x + y*y);
	}
	//dphi/dx
	//dphi/dy
	dboy_dcar[3][3] = 0;//dphi/dz
	
}

void dboydcar_boyknown(double spin, double car[], double boy[], double dboy_dcar[][4]) {//已知两套坐标，快速一点算导数
	double x = car[1], y = car[2], z = car[3];
	double a = spin;
	double r = boy[1], th = boy[2], phi = boy[3];

	dboy_dcar[0][0] = 1;
	dboy_dcar[0][1] = 0;
	dboy_dcar[0][2] = 0;
	dboy_dcar[0][3] = 0;
	dboy_dcar[1][0] = 0;
	dboy_dcar[2][0] = 0;
	dboy_dcar[3][0] = 0;
	dboy_dcar[1][1] = r*r*r*x / (r*r*r*r + a*a*z*z);//dr/dx
	dboy_dcar[1][2] = r*r*r*y / (r*r*r*r + a*a*z*z);//dr/dy
	dboy_dcar[1][3] = r*(a*a + r*r)*z / (r*r*r*r + a*a*z*z);//dr/dz

	if (abs(th - Piby2) < eps) { //th趋于pi/2的极限
		dboy_dcar[2][1] = 0;
		dboy_dcar[2][2] = 0;
		dboy_dcar[2][3] = -1 / r;
	}
	else if (abs(th - 0) < eps) { //th趋于0的极限
		dboy_dcar[2][1] = cos(phi) / sqrt(r*r + a*a);
		dboy_dcar[2][2] = sin(phi) / sqrt(r*r + a*a);
		dboy_dcar[2][3] = 0;
	}
	else {
		dboy_dcar[2][1] = x*cos(th)*cos(th) / (a*a*cos(th)*cos(th)*cos(th)*sin(th) + z*z*tan(th));
		dboy_dcar[2][2] = y*cos(th)*cos(th) / (a*a*cos(th)*cos(th)*cos(th)*sin(th) + z*z*tan(th));
		dboy_dcar[2][3] = -z *sin(th)*sin(th) / (a*a*cos(th)*cos(th)*cos(th)*sin(th) + z*z*tan(th));

	}

	if (abs(x - 0) < eps) {//x趋于0的极限
		dboy_dcar[3][1] = -1 / y;
		dboy_dcar[3][2] = 0;
	}
	else if (abs(y - 0)<eps) {//y趋于0的极限
		dboy_dcar[3][1] = 0;
		dboy_dcar[3][2] = 1 / x;
	}
	else {
		dboy_dcar[3][1] = -y / (x*x + y*y);
		dboy_dcar[3][2] = x / (x*x + y*y);
	}
	//dphi/dx
	//dphi/dy
	dboy_dcar[3][3] = 0;//dphi/dz

}

void highderiv(double spin, double krz_d[], double boy[], double boy_taudot[], double x[][4]) {//坐标对t的s阶导数存在x[s]里
	//写的时候注意保留x[][]为Cartesian坐标下对t的0阶至7阶导数，Gamma为Cartesian坐标下的Gamma

	double dcar_dboy[4][4] = { 0 }, dboy_dcar[4][4] = { 0 };
	double Gamma_boy[4][4][4],Gamma[4][4][4];
	double r = boy[1], th = boy[2];
	dcardboy(spin, boy, x[0], dcar_dboy);//算出cartesian坐标下的0阶导数
	dboydcar_boyknown(spin, x[0], boy, dboy_dcar);
	Christoffel_KRZ(spin, krz_d, r, th, Gamma_boy);
	int al, bt, gm, mu, sg, ro;//指标alpha, beta, gamma, mu, sigma, rho
	
	//转换Christoffel到Cartesian 坐标 {\Gamma^{(Cart)} }^\alpha_{\beta \gamma} = \Gamma^{(Boyer)} ^\mu_{\sigma\rho} \frac{\partial x^{(cart)}^\alpha }{x^{(boyer)} ^\mu } \frac{\partial x^{(boyer)}^\sigma }{x^{(cart)} ^\beta } \frac{\partial x^{(boyer)}^\rho }{x^{(cart)} ^\gamma }
	for (al = 0;al < 4;al++) {
		for (bt = 0;bt < 4;bt++) {
			for (gm = 0;gm < 4;gm++) {
				Gamma[al][bt][gm] = 0;
				for (mu = 0;mu < 4;mu++) {
					for (sg = 0;sg < 4;sg++) {
						for (ro = 0;ro < 4;ro++) {
							Gamma[al][bt][gm] += Gamma_boy[mu][sg][ro] * dcar_dboy[al][mu] * dboy_dcar[sg][bt] * dboy_dcar[ro][gm];
						}
					}
				}
			}
		}
	}

	double car_taudot[4];//速度也换成Cartesian坐标
	for (al = 0;al < 4;al++) {
		car_taudot[al] = 0;
		for (mu = 0;mu < 4;mu++) {
			car_taudot[al] += boy_taudot[mu] * dcar_dboy[al][mu];
		}
	}

	//对t的1阶导数
	x[1][0] = 1;x[1][1] = car_taudot[1] / car_taudot[0]; x[1][2] = car_taudot[2] / car_taudot[0]; x[1][3] = car_taudot[3] / car_taudot[0];

	//对t的2阶导数
	x[2][0] = 0;
	int i;//遍历1到3的指标i
	int nu;
	for (i = 1;i < 4;i++) {
		x[2][i] = 0;
		for (mu = 0;mu < 4;mu++) {
			for (nu = 0;nu < 4;nu++) {
				x[2][i] += (-Gamma[i][mu][nu] * x[1][mu] * x[1][nu] + Gamma[0][mu][nu] * x[1][mu] * x[1][nu] * x[1][i]);
			}
		}
	}

	//对t的3阶导数
	x[3][0] = 0;
	for (i = 1;i < 4;i++) {
		x[3][i] = 0;
		for (mu = 0;mu < 4;mu++) {
			for (nu = 0;nu < 4;nu++) {
				x[3][i] += (-Gamma[i][mu][nu] * (x[2][mu] * x[1][nu] + x[1][mu] * x[2][nu])
					+ Gamma[0][mu][nu] * (x[2][mu] * x[1][nu] * x[1][i] + x[1][mu] * x[2][nu] * x[1][i] + x[1][mu] * x[1][nu] * x[2][i]));
			}
		}
	}

	//对t的4阶导数
	x[4][0] = 0;
	for (i = 1;i < 4;i++) {
		x[4][i] = 0;
		for (mu = 0;mu < 4;mu++) {
			for (nu = 0;nu < 4;nu++) {
				x[4][i] += (-Gamma[i][mu][nu] * (x[3][mu] * x[1][nu] + 2 * x[2][mu] * x[2][nu] + x[1][mu] * x[3][nu])
					+ Gamma[0][mu][nu] * (x[3][mu] * x[1][nu] * x[1][i] + x[1][mu] * x[3][nu] * x[1][i] + x[1][mu] * x[1][nu] * x[3][i]
						+ 2 * x[2][mu] * x[2][nu] * x[1][i] + 2 * x[1][mu] * x[2][nu] * x[2][i] + 2 * x[2][mu] * x[1][nu] * x[2][i]));
			}
		}
	}

	//对t的5阶导数
	x[5][0] = 0;
	for (i = 1;i < 4;i++) {
		x[5][i] = 0;
		for (mu = 0;mu < 4;mu++) {
			for (nu = 0;nu < 4;nu++) {
				x[5][i] += (-Gamma[i][mu][nu] * (x[4][mu] * x[1][nu] + 3 * x[3][mu] * x[2][nu] + 3* x[2][mu]*x[3][nu] + x[1][mu] * x[4][nu])
					+ Gamma[0][mu][nu] * (x[4][mu] * x[1][nu] * x[1][i] + x[1][mu] * x[4][nu] * x[1][i] + x[1][mu] * x[1][nu] * x[4][i]
						+ 3 * x[3][mu] * x[2][nu] * x[1][i] + 3 * x[2][mu] * x[3][nu] * x[1][i] + 3 * x[1][mu] * x[3][nu] * x[2][i] + 3 * x[1][mu] * x[2][nu] * x[3][i] 
						+ 3* x[2][mu]*x[1][nu]*x[3][i] +3* x[3][mu]*x[1][nu]*x[2][i] + 6 * x[2][mu]*x[2][nu]*x[2][i] ));
			}
		}
	}

	//对t的6阶导数
	x[6][0] = 0;
	for (i = 1;i < 4;i++) {
		x[6][i] = 0;
		for (mu = 0;mu < 4;mu++) {
			for (nu = 0;nu < 4;nu++) {
				x[6][i] += (-Gamma[i][mu][nu] * (x[5][mu] * x[1][nu] + 4 * x[4][mu] * x[2][nu] + 6 * x[3][mu]*x[3][nu] + 4 * x[2][mu] * x[4][nu] + x[1][mu] * x[5][nu])
					+ Gamma[0][mu][nu] * (x[5][mu] * x[1][nu] * x[1][i] + x[1][mu] * x[5][nu] * x[1][i] + x[1][mu] * x[1][nu] * x[5][i]
						+ 4 * x[4][mu] * x[2][nu] * x[1][i] + 4 * x[2][mu] * x[4][nu] * x[1][i] + 4 * x[1][mu] * x[4][nu] * x[2][i]
						+ 4 * x[1][mu] * x[2][nu] * x[4][i] + 4 * x[2][mu] * x[1][nu] * x[4][i] + 4 * x[4][mu] * x[1][nu] * x[2][i]
						+ 6 * x[3][mu] * x[3][nu] * x[1][i] + 6 * x[3][mu] * x[1][nu] * x[3][i] + 6 * x[1][mu] * x[3][nu] * x[3][i]
						+ 12* x[2][mu] * x[2][nu] * x[3][i] + 12* x[2][mu] * x[3][nu] * x[2][i] + 12* x[3][mu] * x[2][nu] * x[2][i] ));
			}
		}
	}

	//对t的7阶导数
	x[7][0] = 0;
	for (i = 1;i < 4;i++) {
		x[7][i] = 0;
		for (mu = 0;mu < 4;mu++) {
			for (nu = 0;nu < 4;nu++) {
				x[7][i] += (-Gamma[i][mu][nu] * (x[5][mu] * x[1][nu] + 4 * x[4][mu] * x[2][nu] + 6 * x[3][mu] * x[3][nu] + 4 * x[2][mu] * x[4][nu] + x[1][mu] * x[5][nu])
					+ Gamma[0][mu][nu] * (x[5][mu] * x[1][nu] * x[1][i] + x[1][mu] * x[5][nu] * x[1][i] + x[1][mu] * x[1][nu] * x[5][i]
						+ 4 * x[4][mu] * x[2][nu] * x[1][i] + 4 * x[2][mu] * x[4][nu] * x[1][i] + 4 * x[1][mu] * x[4][nu] * x[2][i]
						+ 4 * x[1][mu] * x[2][nu] * x[4][i] + 4 * x[2][mu] * x[1][nu] * x[4][i] + 4 * x[4][mu] * x[1][nu] * x[2][i]
						+ 6 * x[3][mu] * x[3][nu] * x[1][i] + 6 * x[3][mu] * x[1][nu] * x[3][i] + 6 * x[1][mu] * x[3][nu] * x[3][i]
						+ 12 * x[2][mu] * x[2][nu] * x[3][i] + 12 * x[2][mu] * x[3][nu] * x[2][i] + 12 * x[3][mu] * x[2][nu] * x[2][i]));
			}
		}
	}
}

void radacc(double spin, double x[][4], double acc[]) {//注意，乘上质量比才是真·加速度
	double lev[4][4][4] = { 0 };
	lev[1][2][3] = 1; lev[2][3][1] = 1; lev[3][1][2] = 1;
	lev[3][2][1] = -1; lev[2][1][3] = -1; lev[1][3][2] = -1;//levi-civita记号

	int i, j, k, p, q;//求和的指标
	double I_5dot[4][4] = { 0 }, J_5dot[4][4] = { 0 }, J_6dot[4][4] = { 0 };//0分量放着占位

	//I_5dot[1][1];
	//I_5dot[2][2];
	//I_5dot[3][3];
	//I_5dot[1][2];
	//I_5dot[2][3];
	//I_5dot[1][3];
	for (i = 1;i < 4;i++) {
		for (j = i;j < 4;j++) {
			I_5dot[i][j] = x[5][i] * x[0][j] + 5 * x[4][i] * x[1][j] + 10 * x[3][i] * x[2][j] + 5 * x[2][i] * x[3][j] + x[0][i] * x[5][j];
		}
	}
	I_5dot[2][1] = I_5dot[1][2];
	I_5dot[3][2] = I_5dot[2][3];
	I_5dot[3][1] = I_5dot[1][3];

	//算J的5,6次导数
	int k1, m1, k2, m2;
	for (i = 1;i < 4;i++) {
		for (j = 1;j < 4;j++) {

			if (j == 1) {
				k1 = 2;m1 = 3;
				k2 = 3;m2 = 2;
			}
			else if (j == 2) {
				k1 = 3;m1 = 1;
				k2 = 1;m2 = 3;
			}
			else {
				J_5dot[i][j] = -3.0 / 2.0*spin*x[5][i];
				J_6dot[i][j] = -3.0 / 2.0*spin*x[6][i];//代替delta_j3
				k1 = 1;m1 = 2;
				k2 = 2;m2 = 1;
			}//用来代替levi-civita

			J_5dot[i][j] += (x[5][i] * x[0][k1] * x[1][m1] + x[0][i] * x[5][k1] * x[1][m1] + x[0][i] * x[0][k1] * x[6][m1]
				+ 5 * x[4][i] * x[1][k1] * x[1][m1] + 5 * x[1][i] * x[4][k1] * x[1][m1] + 5 * x[0][i] * x[4][k1] * x[2][m1] + 5 * x[4][i] * x[0][k1] * x[2][m1] + 5 * x[1][i] * x[0][k1] * x[5][m1] + 5 * x[0][i] * x[1][k1] * x[5][m1]
				+ 10 * x[3][i] * x[2][k1] * x[1][m1] + 10 * x[2][i] * x[3][k1] * x[1][m1] + 10 * x[3][i] * x[0][k1] * x[3][m1] + 10 * x[0][i] * x[3][k1] * x[3][m1] + 10 * x[2][i] * x[0][k1] * x[4][m1] + 10 * x[0][i] * x[2][k1] * x[4][m1]
				+ 20 * x[3][i] * x[1][k1] * x[2][m1] + 20 * x[1][i] * x[3][k1] * x[2][m1] + 20 * x[1][i] * x[1][k1] * x[4][m1]
				+ 30 * x[2][i] * x[2][k1] * x[2][m1] + 30 * x[2][i] * x[1][k1] * x[3][m1] + 30 * x[1][i] * x[2][k1] * x[3][m1])
				- (x[5][i] * x[0][k2] * x[1][m2] + x[0][i] * x[5][k2] * x[1][m2] + x[0][i] * x[0][k2] * x[6][m2]
					+ 5 * x[4][i] * x[1][k2] * x[1][m2] + 5 * x[1][i] * x[4][k2] * x[1][m2] + 5 * x[0][i] * x[4][k2] * x[2][m2] + 5 * x[4][i] * x[0][k2] * x[2][m2] + 5 * x[1][i] * x[0][k2] * x[5][m2] + 5 * x[0][i] * x[1][k2] * x[5][m2]
					+ 10 * x[3][i] * x[2][k2] * x[1][m2] + 10 * x[2][i] * x[3][k2] * x[1][m2] + 10 * x[3][i] * x[0][k2] * x[3][m2] + 10 * x[0][i] * x[3][k2] * x[3][m2] + 10 * x[2][i] * x[0][k2] * x[4][m2] + 10 * x[0][i] * x[2][k2] * x[4][m2]
					+ 20 * x[3][i] * x[1][k2] * x[2][m2] + 20 * x[1][i] * x[3][k2] * x[2][m2] + 20 * x[1][i] * x[1][k2] * x[4][m2]
					+ 30 * x[2][i] * x[2][k2] * x[2][m2] + 30 * x[2][i] * x[1][k2] * x[3][m2] + 30 * x[1][i] * x[2][k2] * x[3][m2]);
			J_6dot[i][j]+= (x[6][i] * x[0][k1] * x[1][m1] + x[0][i] * x[6][k1] * x[1][m1] + x[0][i] * x[0][k1] * x[7][m1]
				+ 6 * x[6][i] * x[1][k1] * x[1][m1] + 6 * x[1][i] * x[6][k1] * x[1][m1] + 6 * x[0][i] * x[5][k1] * x[2][m1] + 6 * x[5][i] * x[0][k1] * x[2][m1] + 6 * x[1][i] * x[0][k1] * x[6][m1] + 6 * x[0][i] * x[1][k1] * x[6][m1]
				+ 15 * x[4][i] * x[2][k1] * x[1][m1] + 15 * x[2][i] * x[4][k1] * x[1][m1] + 15 * x[4][i] * x[0][k1] * x[3][m1] + 15 * x[0][i] * x[4][k1] * x[3][m1] + 15 * x[2][i] * x[0][k1] * x[5][m1] + 15 * x[0][i] * x[2][k1] * x[5][m1]
				+ 20 * x[3][i] * x[3][k1] * x[1][m1] + 20 * x[3][i] * x[0][k1] * x[4][m1] + 20 * x[1][i] * x[3][k1] * x[4][m1]
				+ 30 * x[4][i] * x[1][k1] * x[2][m1] + 30 * x[1][i] * x[4][k1] * x[2][m1] + 30 * x[1][i] * x[1][k1] * x[5][m1]
				+ 60 * x[3][i] * x[2][k1] * x[2][m1] + 60 * x[2][i] * x[3][k1] * x[2][m1] + 60 * x[3][i] * x[1][k1] * x[3][m1]
				+ 60 * x[1][i] * x[3][k1] * x[3][m1] + 60 * x[2][i] * x[1][k1] * x[4][m1] + 60 * x[1][i] * x[2][k1] * x[4][m1]
				+ 90 * x[2][i] * x[2][k1] * x[3][m1])
				- (x[6][i] * x[0][k2] * x[1][m2] + x[0][i] * x[6][k2] * x[1][m2] + x[0][i] * x[0][k2] * x[7][m2]
					+ 6 * x[6][i] * x[1][k2] * x[1][m2] + 6 * x[1][i] * x[6][k2] * x[1][m2] + 6 * x[0][i] * x[5][k2] * x[2][m2] + 6 * x[5][i] * x[0][k2] * x[2][m2] + 6 * x[1][i] * x[0][k2] * x[6][m2] + 6 * x[0][i] * x[1][k2] * x[6][m2]
					+ 15 * x[4][i] * x[2][k2] * x[1][m2] + 15 * x[2][i] * x[4][k2] * x[1][m2] + 15 * x[4][i] * x[0][k2] * x[3][m2] + 15 * x[0][i] * x[4][k2] * x[3][m2] + 15 * x[2][i] * x[0][k2] * x[5][m2] + 15 * x[0][i] * x[2][k2] * x[5][m2]
					+ 20 * x[3][i] * x[3][k2] * x[1][m2] + 20 * x[3][i] * x[0][k2] * x[4][m2] + 20 * x[1][i] * x[3][k2] * x[4][m2]
					+ 30 * x[4][i] * x[1][k2] * x[2][m2] + 30 * x[1][i] * x[4][k2] * x[2][m2] + 30 * x[1][i] * x[1][k2] * x[5][m2]
					+ 60 * x[3][i] * x[2][k2] * x[2][m2] + 60 * x[2][i] * x[3][k2] * x[2][m2] + 60 * x[3][i] * x[1][k2] * x[3][m2]
					+ 60 * x[1][i] * x[3][k2] * x[3][m2] + 60 * x[2][i] * x[1][k2] * x[4][m2] + 60 * x[1][i] * x[2][k2] * x[4][m2]
					+ 90 * x[2][i] * x[2][k2] * x[3][m2]);
		}
	}

	//算加速度
	acc[0] = 0;
	for (j = 1;j < 4;j++) {
		acc[j] = 8.0 / 15.0*spin*J_5dot[3][j];
		for (k = 1;k < 4;k++) {
			acc[j] += -2.0 / 5.0*I_5dot[j][k] * x[0][k];
			for (p = 1;p < 4;p++) {
				for (q = 1;q < 4;q++) {
					acc[j] += 16.0 / 45.0 * lev[j][p][q] * J_6dot[p][k] * x[0][q] * x[0][k] + 32.0 / 45.0 * lev[j][p][q] * J_5dot[p][k] * x[0][k] * x[1][q]
						+16.0/45.0 * lev[p][q][j] * J_5dot[k][p] * x[0][q] *x[1][k] - 16.0/45.0 * lev[p][q][k] * J_5dot[j][p] * x[0][q] * x[1][k];
				}
			}
		}
	}
}
//一个以前用的RK4算法
//for (i = 0; i < 1000000;i++) {
//	/*if (r > 30) {
//	system("pause");
//	}*/
//	if (var[1] < horizon + 0.01) {
//		printf("Fell into horizon\n");
//		break;
//	}
//	for (;;) {
//		if (var[1] < horizon + 0.01) {
//			printf("Fell into horizon\n");
//			break;
//		}
//		t = var[0];	r = var[1];	th = var[2];	phi = var[3];	ur = var[4];	uth = var[5];	ut = var[6];	uphi = var[7];
//
//		/******************** CALCULATE k1 *****************/
//		metric_KRZ(spin, krz_d, r, th, g);
//		Christoffel_KRZ(spin, krz_d, r, th, Gamma);
//		k1[0] = ut;	varnew1[0] = var[0] + 0.5*dtau*k1[0];
//		k1[3] = uphi;	varnew1[3] = var[3] + 0.5*dtau*k1[3];
//		k1[2] = uth;	varnew1[2] = var[2] + 0.5*dtau*k1[2];
//		k1[1] = ur;		varnew1[1] = var[1] + 0.5*dtau*k1[1];
//
//
//
//		u[0] = ut;
//		u[1] = ur;
//		u[2] = uth;
//		u[3] = uphi;
//
//
//
//
//		F_r = 0;
//		F_theta = 0;
//		F_t = 0;
//		F_phi = 0;
//		ita = 0;
//		for (ii = 0;ii < 4;ii++) {
//			for (jj = 0;jj < 4;jj++) {
//
//				F_t -= 0.5*Gamma[0][ii][jj] * u[ii] * u[jj];
//				F_r -= 0.5*Gamma[1][ii][jj] * u[ii] * u[jj];
//				F_theta -= 0.5*Gamma[2][ii][jj] * u[ii] * u[jj];
//				F_phi -= 0.5*Gamma[3][ii][jj] * u[ii] * u[jj];
//				ita += g[ii][jj] * u[ii] * u[jj];
//			}
//		}
//		//if (ita > -0.8) {
//		//system("pause");
//		//break;
//		//ita = -0.8;
//		//}
//		E = -g[0][0] * ut - g[0][3] * uphi;
//		Lz = g[0][3] * ut + g[3][3] * uphi;
//		Q = uth*g[2][2] * uth*g[2][2] + cos(th)*cos(th)*(spin*spin*(mu*mu - E*E) + Lz*Lz / sin(th) / sin(th));
//
//		k1[4] = F_r;	varnew1[4] = var[4] + 0.5*dtau*k1[4];
//		k1[5] = F_theta;	varnew1[5] = var[5] + 0.5*dtau*k1[5];
//		k1[6] = F_t;	varnew1[6] = var[6] + 0.5*dtau*k1[6];
//		k1[7] = F_phi;	varnew1[7] = var[7] + 0.5*dtau*k1[7];
//		/*********************finished calculate k1********************/
//
//		/******************** CALCULATE k2 *****************/
//		metric_KRZ(spin, krz_d, varnew1[1], varnew1[2], g);
//		Christoffel_KRZ(spin, krz_d, varnew1[1], varnew1[2], Gamma);
//
//		ut = varnew1[6];	k2[0] = varnew1[6];	varnew2[0] = var[0] + 0.5*dtau*k2[0];
//		ur = varnew1[4];	k2[1] = varnew1[4];	varnew2[1] = var[1] + 0.5*dtau*k2[1];
//		uth = varnew1[5];	k2[2] = varnew1[5];	varnew2[2] = var[2] + 0.5*dtau*k2[2];
//		uphi = varnew1[7];	k2[3] = varnew1[7];	varnew2[3] = var[3] + 0.5*dtau*k2[3];
//
//		u[0] = ut;	u[1] = ur;	u[2] = uth;	u[3] = uphi;
//
//		F_r = 0;
//		F_theta = 0;
//		F_t = 0;
//		F_phi = 0;
//		for (ii = 0;ii < 4;ii++) {
//			for (jj = 0;jj < 4;jj++) {
//
//				F_t -= 0.5*Gamma[0][ii][jj] * u[ii] * u[jj];
//				F_r -= 0.5*Gamma[1][ii][jj] * u[ii] * u[jj];
//				F_theta -= 0.5*Gamma[2][ii][jj] * u[ii] * u[jj];
//				F_phi -= 0.5*Gamma[3][ii][jj] * u[ii] * u[jj];
//			}
//		}
//		k2[4] = F_r;	varnew2[4] = var[4] + 0.5*dtau*k2[4];
//		k2[5] = F_theta;	varnew2[5] = var[5] + 0.5*dtau*k2[5];
//		k2[6] = F_t;	varnew2[6] = var[6] + 0.5*dtau*k2[6];
//		k2[7] = F_phi;	varnew2[7] = var[7] + 0.5*dtau*k2[7];
//
//		/*********************finished calculate k2********************/
//
//
//		/*********************CALCULATE k3********************/
//		metric_KRZ(spin, krz_d, varnew2[1], varnew2[2], g);
//		Christoffel_KRZ(spin, krz_d, varnew2[1], varnew2[2], Gamma);
//
//		ut = varnew2[6];	k3[0] = varnew2[6];	varnew3[0] = var[0] + dtau*k3[0];
//		ur = varnew2[4];	k3[1] = varnew2[4];	varnew3[1] = var[1] + dtau*k3[1];
//		uth = varnew2[5];	k3[2] = varnew2[5];	varnew3[2] = var[2] + dtau*k3[2];
//		uphi = varnew2[7];	k3[3] = varnew2[7];	varnew3[3] = var[3] + dtau*k3[3];
//
//		u[0] = ut;	u[1] = ur;	u[2] = uth;	u[3] = uphi;
//
//		F_r = 0;
//		F_theta = 0;
//		F_t = 0;
//		F_phi = 0;
//		for (ii = 0;ii < 4;ii++) {
//			for (jj = 0;jj < 4;jj++) {
//
//				F_t -= 0.5*Gamma[0][ii][jj] * u[ii] * u[jj];
//				F_r -= 0.5*Gamma[1][ii][jj] * u[ii] * u[jj];
//				F_theta -= 0.5*Gamma[2][ii][jj] * u[ii] * u[jj];
//				F_phi -= 0.5*Gamma[3][ii][jj] * u[ii] * u[jj];
//			}
//		}
//		k3[4] = F_r;	varnew3[4] = var[4] + dtau*k3[4];
//		k3[5] = F_theta;	varnew3[5] = var[5] + dtau*k3[5];
//		k3[6] = F_t;	varnew3[6] = var[6] + dtau*k3[6];
//		k3[7] = F_phi;	varnew3[7] = var[7] + dtau*k3[7];
//
//		/*********************finished calculate k3********************/
//
//
//		/*********************CALCULATE k4********************/
//		metric_KRZ(spin, krz_d, varnew3[1], varnew3[2], g);
//		Christoffel_KRZ(spin, krz_d, varnew3[1], varnew3[2], Gamma);
//
//		ut = varnew3[6];	k4[0] = varnew3[6];
//		ur = varnew3[4];	k4[1] = varnew3[4];
//		uth = varnew3[5];	k4[2] = varnew3[5];
//		uphi = varnew3[7];	k4[3] = varnew3[7];
//
//		u[0] = ut;	u[1] = ur;	u[2] = uth;	u[3] = uphi;
//
//		F_r = 0;
//		F_theta = 0;
//		F_t = 0;
//		F_phi = 0;
//		for (ii = 0;ii < 4;ii++) {
//			for (jj = 0;jj < 4;jj++) {
//
//				F_t -= 0.5*Gamma[0][ii][jj] * u[ii] * u[jj];
//				F_r -= 0.5*Gamma[1][ii][jj] * u[ii] * u[jj];
//				F_theta -= 0.5*Gamma[2][ii][jj] * u[ii] * u[jj];
//				F_phi -= 0.5*Gamma[3][ii][jj] * u[ii] * u[jj];
//			}
//		}
//		k4[4] = F_r;
//		k4[5] = F_theta;
//		k4[6] = F_t;
//		k4[7] = F_phi;
//		/****************************finished calculate k4********************/
//
//		if (fabs(fabs(dtau*(k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0]) / 6.0) - dt) < tol) break;
//		else {
//			dtau = fabs(6.0*dt / (k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0]));
//		}
//	}
//
