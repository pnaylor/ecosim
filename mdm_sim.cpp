// ************************************************************************************
//        Mussel Disturbance Model (MDM) based on the Forest Fire Model (FFM)
//        Lattice simulation with cluster algorithm
//        written by Frederic Guichard, Princeton University, last update 02/2003
//
//	  Adapted as model for multispecies assembly 
//	  Multi species assembly model written by Pradeep Pillai, McGill Unviversity
//	  Last update 09/2005
//
//        Requires cokus.h header file in the include path for the cokus random number generator
// On Linux, compile with:
// g++ -lm -O3 -g mdm_sim.cc -o essai `pkg-config --cflags gtk+-2.0` `pkg-config --libs gtk+-2.0`
// the command calls the pkg-config programs which takes care of including the required libraries
// ************************************************************************************

#include <fstream>
#include <iostream>
using namespace std;
#include <omp.h>
#include "cokus.h"
//#include "RandMT.h"  // random number generator (not in standard C distribution)
//using namespace RandMT;

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vector>


#define T 300                        //(originally 50000) How many time steps to output (T>=1)
#define BEGIN 100                         // when to start outputs
#define MIGRANTPOOL 1000    //Number of new species in pool(#potential migrants)
int step;

int Simulation_time = 0;

#include "mdm_model.h"   // code containing the actual simulation. This file would be independent with the addition of a main function and of the global variables declared above.





int main(int argc, char *argv[])
{
	cout << omp_get_num_procs() << " processors available" << endl;

  FILE *fptr;
  char filename[]="richness.txt";
  fptr=fopen(filename, "w");
  fclose(fptr);
  
  init_parameters();
  init_matcover();
  seedrand(); // initialize the cokus.h (required file) random number generator with the clock
  
  initlattice();
  step=0;
  ab2=0;
  
	//cout << NSITES << endl;


  for(int pradeep=0; pradeep<T+BEGIN; pradeep++)
  {
     update_lattice(step);	// run the simulation 1 time step (NSITES cell updates)
     
     //call functions if the simulation has completed BEGIN time steps
     if (step>=BEGIN && step<T+BEGIN)
     {
      //cluster();  // will run the cluster algorithm
      //outdens(step); // will append data to the time series output file
      //ab2cum=ab2cum+ab2tot; //increment the total mussel cover value
 

       migrant_abundance();  // IMPORTANT! required every time step because several 
                             // functions rely on this to function
   
       //update_spDuration(); //update the spDuration[] array


    
       int rich=richness();
       //printf("Species Richness at time step %d is: %d\n", step, rich);
	   //cout << "Species Richness at time step " << step << " is: " << rich << endl;
	   sumRichness+=rich;
       cycleCount++;

    
//       fout=fopen(filename, "a");
//       fprintf(fout, "%d \t %d \n",step,rich);
//       fclose(fout);

     }
  
     step++;
	 cout << ' ' << step;

  }


	//cout << "finished";
//--------------------------------------------------------------------------------


  
  
  //out(); // write output to file
  //outmap(); // output a snapshot of the final lattice
  


  //------call fxns to output abundance for species and ranks abundance----

  //DELETED ALL FUNCTION CALLS IN THIS STRIPPED-DOWN VERSION OF CODE

  //-------end output to files------

  
  return(0);  
}















