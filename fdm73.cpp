#include<stdio.h>
#include<iostream>
#include<math.h>
#include<fstream>
#include <algorithm>
using namespace std;
int main()
{
    int N=61; // number of points
	float xi[200], ai[200]; // geom
	float rho[200], tem[200], vel[200], pre[200]; // flow variables
	float rhom[200], fm[200], rm[200], um[200], tm[200]; // 
	float delt[200], drdt[200], dudt[200], dtdt[200], drdtm[200], dudtm[200];
	float dtdtm[200], drdtav[200], dudtav[200], dtdtav[200];
	float rt[1000], tt[1000], pt[1000], amt[1000];
	float drar[1000], drau[1000], drat[1000];
	float ma[200], me[200];
	float dt=0.0;

	float gamma = 1.4;
	float cfl = 0.5;
	int i, j, k;
	//ofstream OutFile;
	for ( i = 1; i <= N; i++)
	{
		//define geometry
		xi[i] = (i - 1) / ((N - 1) / 3.0);
		ai[i] = 1.0 + 2.2*(xi[i] - 1.5)*(xi[i] - 1.5);
		//init solution
		rho[i] = 1. - 0.3146*xi[i];
		tem[i] = 1. - 0.2314*xi[i];
		vel[i] = (0.1 + 1.09*xi[i])*(sqrt(tem[i]));
		pre[i] = rho[i] * tem[i];
		ma[i] = vel[i] / sqrt(tem[i]);
		me[i] = rho[i] * vel[i] * ai[i];

		//printf("r=%f,a=%f\n",rho[i],ai[i]);
	}
	float deltx = 3.0 / (N - 1);
	float deltt = 0.5;
	//loop

	for (j = 1; j < 2000; j++)
	{
		//find a suitable timestep
		for (i = 1; i <= N; i++)
		{
			delt[i] = cfl*deltx / (sqrt(tem[i]) + vel[i]);
		}
		 dt = delt[2];
		for (i = 3; i < N; i++)
		{
			dt = min(dt, delt[i]);
		}
		//fdm
		for (k = 2; k < N; k++)
		{
			drdt[i] = (-rho[i] * (vel[i + 1] - vel[i]) - rho[i] * vel[i] * (log(ai[i + 1]) - \
				log(ai[i])) - vel[i] * (rho[i + 1] - rho[i])) / deltx;

			dudt[i] = (-vel[i] * (vel[i + 1] - vel[i]) - 1.0 / gamma*(tem[i + 1] - tem[i] + tem[i] / rho[i] * \
				(rho[i + 1] - rho[i]))) / deltx;

			dtdt[i] = (-vel[i] * (tem[i + 1] - tem[i]) - (gamma - 1.0)*tem[i] * (vel[i + 1] - \
				vel[i] + vel[i] * (log(ai[i + 1]) - log(ai[i])))) / deltx;

			//at inlet
			rhom[1] = 1.0;
			um[1] = 2.*um[2] - um[3];
			tm[1] = 1.0;

			//all points
			rhom[i] = rho[i] + drdt[i] * deltt;
			um[i] = vel[i] + dudt[i] * deltt;
			tm[i] = tem[i] + dtdt[i] * deltt;
			
			drdtm[i] = (-rhom[i] * (um[i] - um[i - 1]) - rhom[i] * um[i] * (log(ai[i]) - \
				log(ai[i - 1])) - um[i] * (rhom[i] - rhom[i - 1])) / deltx;
			dudtm[i] = (-um[i] * (um[i] - um[i - 1]) - 1. / gamma*(tm[i] - tm[i - 1] + \
				tm[i] / rhom[i] * (rhom[i] - rhom[i - 1]))) / deltx;
			dtdtm[i] = (-um[i] * (tm[i] - tm[i - 1]) - (gamma - 1.)*tm[i] * (um[i] - um[i - 1] + \
				um[i] * (log(ai[i]) - log(ai[i - 1])))) / deltx;
			//average d()/dt
			drdtav[i] = (drdt[i] + drdtm[i]) / 2.;
			dudtav[i] = (dudt[i] + dudtm[i]) / 2.;
			dtdtav[i] = (dtdt[i] + dtdtm[i]) / 2.;
			//new value
			rho[i] = rho[i] + drdtav[i] * deltt;
			vel[i] = vel[i] + dudtav[i] * deltt;
			tem[i] = tem[i] + dtdtav[i] * deltt;

			pre[i] = rho[i] * tem[i];
			ma[i] = vel[i] / sqrt(tem[i]);
			me[i] = rho[i] * vel[i] * ai[i];

		}
		//at boundary
            //inl
		vel[1] = 2.*vel[2] - vel[3];
		rho[1] = 1.;
		tem[1] = 1.;
		pre[1] = rho[1] * tem[1];
		ma[1] = vel[1] / sqrt(tem[1]);
		me[1] = rho[1] * vel[1] * ai[1];
			//out
		vel[N] = 2.*vel[N - 1] - vel[N - 2];
		rho[N] = 2.*rho[N - 1] - rho[N - 2];
		tem[N] = 2.*tem[N - 1] - tem[N - 2];
		pre[N] = rho[N] * tem[N];
		ma[N]  = vel[N] / sqrt(tem[N]);
		me[N]  = rho[N] * vel[N] * ai[N];

	}  // end loop

	//output to files
	ofstream outfile;
	outfile.open("flow.plt");
	outfile <<"T= \"073\" " << endl;
	outfile << "variables = \"x\",\"a\",\"rho \",\"u \",\"p \",\"ma\"" << endl;
	outfile << "zone f=point" << endl;
	for (i = 1; i <= N; i++)
	{
		outfile << xi[i]<< "	"
			    << ai[i]<< "	"
			    <<rho[i]<< "	"
				<<vel[i]<< "	"
				<<pre[i]<< "	"
				<< ma[i]<< endl;
	}
}
