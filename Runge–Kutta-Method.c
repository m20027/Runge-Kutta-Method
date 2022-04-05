/*ルンゲ・クッタ法のプログラム例*/
#include<stdio.h>
#include<math.h>

#define M_STEP 10001 /*ステップ数の上限*/
#define N_var 3 /*変数の数*/
#define T_s 0.0
#define T_e 100.0
#define Gamma 0.20833333333
#define R0 2.5
#define N 127443493.0


void func(double t, double x[], double  f[])/*f(t,x)を定義する*/
{
	double lambda;
	lambda=Gamma*R0/N*x[1];
	f[0]=-lambda*x[0];
	f[1]=lambda*x[0]-Gamma*x[1];
	f[2]=Gamma*x[1];
}

void initialize(double x[])/*初期値を与える*/
{
	x[0]=N-10.0;
	x[1]=10.0;
	x[2]=0.0;
}

int main()
{
	int n,i,n_step;
	double h,t[M_STEP],x[M_STEP][N_var];/*解の近似値*/
	double fk1[N_var],fk2[N_var],fk3[N_var],fk4[N_var],y[N_var],liap,liap0;/*中間変数*/

	initialize(x[0]);/*初期条件の設定*/
	n_step=1000;
	h=(T_e-T_s)/n_step;

	for(n=0;n<n_step;n++)/*ルンゲ・クッタ法の計算*/
	{
		t[n]=T_s+n*h;
		
		func(t[n],x[n],fk1);
		
		for(i=0;i<N_var;i++)
		{
			y[i]=x[n][i]+0.5*h*fk1[i];
		}
		
		func(t[n]+0.5*h,y,fk2);

		for(i=0;i<N_var;i++)
		{
			y[i]=x[n][i]+0.5*h*fk2[i];
		}

		func(t[n]+0.5*h,y,fk3);
		
		for(i=0;i<N_var;i++)
		{
			y[i]=x[n][i]+h*fk3[i];
		}

		func(t[n]+h,y,fk4);

		for(i=0;i<N_var;i++)
		{
			x[n+1][i]=x[n][i]+h*(fk1[i]+2.0*fk2[i]+2.0*fk3[i]+fk4[i])/6.0;
		}
	}
	for(i=0;i<N_var;i++)/*x(T_e)を出力*/
	{
		printf("x_%d(%lf)=%lf\n",i,T_e,x[n_step][i]);
	}
	return 0;
}
