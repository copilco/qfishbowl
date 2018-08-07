/*
 *  testLaser.cpp
 *  
 *
 *  Created by Camilo Ruiz Mendez on 28/04/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "laser.h"

int main()
{
	FILE *out0;
	
	out0=fopen("testlaserPi.txt","w");
	laser laserx(0.05, 1.e14, 0.057, 10, pi/2.);

	
	for(int ktime=0;ktime<laserx.nmaxt;ktime++)
	{
		fprintf(out0,"%e %e %e %e \n",laserx.e.t[ktime],laserx.e.f[ktime][0] ,laserx.av.f[ktime][0]);
	}
	
}

