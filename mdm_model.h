
// ************************************************************************************
//        Mussel Disturbance Model (MDM) based on the Forest Fire Model (FFM)
//        Lattice simulation with cluster algorithm
//        written by Frederic Guichard, Princeton University, last update 02/2002
//
//
//	  Adapted as model for multispecies assembly 
//	  Multi species assembly model written by Pradeep Pillai, McGill Unviversity
//	  Last update 09/2005
//
//        Requires cokus.h header file for the cokus random number generator
// ************************************************************************************

//#include "RandMT.h"

// *********** Global variables (not recommended but faster) ********//


int matstate[NSITES];                                                               // main lattice state matrix
int S,Sl,apos,psize,i,j,k,n,loop_size;
int newNeighbor,xNgh,yNgh,Xmoins,Xplus,Ymoins,Yplus;                            // for neighborhood function
//int simnbr1,simnbr2;                                                                // simulation loop variables



float d,pnew,diff0,diff2,alpha4ini,alpha0,alpha4;                    // parameters
int max_dist_size,delta0,delta2;
double ab2,ab3,ab2cum,condprob,ab2tot,ab2_array[T],t_array[T];                                             // output variables
//int l;
tm *ptm;                                                                           // clock variable (depends on time.h)
time_t *cur_time;                                                                  // clock variable (depends on time.h)
//************New variables for extended model*********


// *********************** End *************************************

// ******** Function declaration **************
void out();                                                                  // write spatial analysis results to a file
void outmap();                                                           // write the lattice to a file
void outdens(int);                                                         // write the time series to a file
void cluster();                                                           // cluster algorithm
int XYNeighbor8(int,int);                                            // neighborhood function
void update(int);                                                       // apply transition rules on a cell
void update_gradient(int);                                                       // apply transition rules on a cell with a gradient
void update_2states_gradient(int);                                                       // apply transition rules on a cell with a gradient
void update_2states(int);                                                       // apply transition rules with 2 state
void seedrand();                                                       // seed the random number generator
void setoutfiles(int,int);                                             // initialize output files
void  initlattice();                                                        // set inital conditions on the lattice
void new_disturbance(int);
void new_disturbance_size(int,int);
void new_colonization(int);
int x_from_pos(int);
int y_from_pos(int);
int pos_from_xy(int,int);
int periodic_pos(int);
void update_lattice(int);
void init_parameters();
void input();
// ******** End ************//


//************New variables*****************
bool graphics;
bool autostart;
bool parTest;
int iGr;
int iAs;
int par;
int model=0;
float paramValue;
float gamma0;
float gamma2;
int migrantAbund[MIGRANTPOOL];  //abundance of migrants
/*
int abundanceArray[MIGRANTPOOL];  //array of species abundance
int spDuration[T][MIGRANTPOOL];    //duration time for species
int Duration[T*MIGRANTPOOL];  //duration time for sp. ranked from low to high
int extinctionSize[T];
*/

int sumRichness=0;
int cycleCount=0;
char filename[]="richness.txt";

/*
char filename2[]="Species_Abundance.txt", filename3[]="RankAbundance.txt";
char filename4[]="Richness-Plot.txt";
char filename5[]="SpeciesDuration.txt";
FILE *fout;
FILE *fptr2, *fptr3, *fptr4, *fptr5;
*/

int runNumber;

/*
// ----turnover varibles------
FILE *fptr9a, *fout9a;
char filename9a[]="turnover-plot.txt";
*/

int matcover[NSITES];

/*
float totalTurnover=0;   // sum total of the mussel cover turnover rates from each time step
int turnoverCount=0;
int turnoverTotal=0;
*/
// -----------------------------



//************New fxns for extended model************************
void newSpecies_colonization(int aRate);
int XYNeighbour8_newSpecies(int aPos);
int neighborSpecies(int aPos);
int select_Neighbor(int);
int mussel_abundance();
void migrant_abundance();
//void rank_abundance();  //rank species abundance from most to least abundance
int richness(); //measures species richness

/*
void update_spDuration();  //updates species duration
void init_spDuration(); //initializes spDuration array
int speciesCount(int);   //returns number of particular species
void out_abundance();
void out_rankAbund();
void out_richnessPlot();
void out_Duration();
void out_extinction();
void out_musselCover();
void out_richnessHabitat();
void init_matcover();
void musselTurnOverUpdate();
void out_aveMusselTurnover();
void sp_area_curve();
*/

//RandMT r;


void init_parameters(){
  /*ifstream input("mdm_par.txt");
  if(!input)
  {
      cout<<"\n\n\tCouldn't open file\n";
  }*/
  //input>>delta0>>alpha0>>delta2>>alpha4>>gamma0>>gamma2>>iGr>>iAs>>par>>model>>runNumber;
       
  
  graphics = iGr;
  autostart = iAs;
  parTest = par;

  paramValue=alpha4;  //default for parValue
  if (autostart)
  {
	if(parTest)
	{	
	   /*alpha0=alpha0/20;*/
           paramValue=alpha0;
	}
	else
	{
	   /*alpha4=alpha4/20;*/
	   paramValue=alpha4;
	}
  }

  
  //delta0=1; // number of new disturbances / time step
  //delta2=10; // number of global colonization attempts / time step
  //alpha0=1; // probability of disturbance spreading
  //alpha4ini=0.6; // initial probability of recovery spreading
  //gamma0=1;  //number of new species per time step
  //gamm2=0.1; //probability of new species spreading
  //alpha4=alpha4ini; // probability of recovery spreading
  max_dist_size=1;//(int) (SIZE/20);
  loop_size=1;
 
}


void update_lattice(int tstep){
  int t;
  //cout << NSITES << endl;
  for(t=0;t<NSITES;t++){
	  //cout << t << ' ';
	int  c = -1;
	while(c < 0 || c >= NSITES)
	{
		c =(int) (randomMT()*(NSITES-1)); // select a cell at random
		//cout << c << ' ';
	}
    update(c);
  }
  //cout << endl;
  
  new_disturbance_size(delta0,max_dist_size); //potentially create a new disturbance
  new_colonization(delta2); //potentially add one or a small number of mussels to a few random sites
  
  //determine global colonization rate for a new sp.

  int colonznNumber=0; //number of cells to be 
                       //colonized by a new species
  int extra=0;

  if(gamma0<=1)
  {
   	  if (randomMT()<=gamma0)
  		  extra=1;
  	  else
  		  extra=0;
  } 
  else
  {
	  if(randomMT()< (gamma0-(int)gamma0) )
		  extra=1;
	  else
		  extra=0;
  }
  colonznNumber= extra + (int)gamma0;

  newSpecies_colonization(colonznNumber);

//    newSpecies_colonization(gamma0);

}           



void seedrand(){
  cur_time = new time_t;
  ptm = new tm;
  time(cur_time);
  ptm = localtime(cur_time);
  int time=int(ptm->tm_hour)*int(ptm->tm_min)*int(ptm->tm_sec);
  while(int (time % 2) == 0) time++;
  seedMT(time);
  }





void  initlattice(){
  init_parameters();
      for (i=0;i<NSITES;++i){
        matstate[i]=2;

      }
      ab2cum=0;

      if(delta0==0){   // initial disturbance front if delta0=0
	for (i=0;i<SIZE;i++){
	  matstate[i]=3;
	}
      }

}




void new_colonization(int aRate){
  for(j=0;j<aRate;j++){
    int aPos=(int) (SIZE*SIZE*randomMT());
    if(matstate[aPos]!=3) matstate[aPos]=2;
  }
}


void newSpecies_colonization(int aRate)
{
	for (j=0;j<aRate;j++)
	{
		int aPos= (int)(SIZE*SIZE*randomMT());
                int newspec=(int)(MIGRANTPOOL*randomMT())+100; //speciesNumber
		if (matstate[aPos]==0 || matstate[aPos]>newspec)
	             matstate[aPos]= newspec;                
			
	}
}



void new_disturbance(int aRate){
  for(j=0;j<aRate;j++){
    int aPos=(int) (SIZE*SIZE*randomMT());
    if(matstate[aPos]==2) matstate[aPos]=3;
  }
}

void new_disturbance_size(int aRate, int max_size){
  for(int k=0;k<aRate;k++){
    int dist_size=1;  //(int) (max_size*randomMT());
    int aPos=(int) (SIZE*SIZE*randomMT());
    int thePos;

    for (i=0;i<dist_size;i++){
      for (j=0;j<dist_size;j++){
	thePos=pos_from_xy(periodic_pos(x_from_pos(aPos)+i),periodic_pos(y_from_pos(aPos)+j));
	if(matstate[thePos]==2) matstate[thePos]=3;
      }
    }
  }
}


// update() function displays hierarchical rules for spread of new species.
// Trade-off between higher fecundity(prob of spreading to neighbour cell) and lower
// competitive ability(can be replaced by slower spreading species)

void update(int aPos){      
  //cout << aPos << endl;
  // ************ Transition rules ****************//
  if(matstate[aPos]==2 && XYNeighbor8(aPos,3)>0 && randomMT()<=alpha0) matstate[aPos]=3;
  else if(matstate[aPos]==3) matstate[aPos]=0;
  else if (matstate[aPos]!=3)
        {  
           int ngh=select_Neighbor(aPos); 
           if(ngh==2 && alpha4>=randomMT())
              matstate[aPos]=2;
           else 
	      if(ngh>99 && (matstate[aPos]==0 || matstate[aPos]>ngh) && ((ngh-99)*((float)1/MIGRANTPOOL))>=randomMT()) //trade off between fecundity and competitiveness
	          matstate[aPos]=ngh;
         }
  
  // ************** end ****************//

}



void update_gradient(int aPos){
  
  // ************ Transition rules ****************//
  if(matstate[aPos]==2 && XYNeighbor8(aPos,3)>0 && randomMT()<=((double)(aPos%SIZE)/(double)(SIZE))) matstate[aPos]=3;
  else if(matstate[aPos]==0 && (0.125*XYNeighbor8(aPos,2)*alpha4)>=randomMT()) matstate[aPos]=2;
  else if(matstate[aPos]==3) matstate[aPos]=0;
  // ************** end ****************//

}

void update_2states_gradient(int aPos){
  
  // ************ Transition rules ****************//
  if(matstate[aPos]==2 && XYNeighbor8(aPos,3)>0 && randomMT()<=((double)(aPos%SIZE)/(double)(SIZE))) matstate[aPos]=3;
  else if(matstate[aPos]==3 && (0.125*XYNeighbor8(aPos,2)*alpha4)>=randomMT()) matstate[aPos]=2;
  // ************** end ****************//

}

void update_2states(int aPos){
  
  // ************ Transition rules ****************//
  if(matstate[aPos]==2 && XYNeighbor8(aPos,3)>0 && randomMT()<=alpha0) matstate[aPos]=3;
  else if(matstate[aPos]==3 && (0.125*XYNeighbor8(aPos,2)*alpha4)>=randomMT()) matstate[aPos]=2;
  // ************** end ****************//

}


int XYNeighbor8(int aPos, int value)
{

  xNgh =aPos%SIZE;
  yNgh =(int) (aPos>>DIM); 

  // Define neighborhood coordinates
  Xmoins=xNgh-1;
  Xplus=xNgh+1;
  Ymoins=yNgh-1;
  Yplus=yNgh+1;

  // Apply periodic boundary conditions
  if((xNgh-1)<0) Xmoins=SIZE-1;
  else if((xNgh+1)>(SIZE-1)) Xplus=0;
  if((yNgh-1)<0) Ymoins=SIZE-1;
  else if((yNgh+1)>(SIZE-1)) Yplus=0;

  // Compute and return the number of neighbors in state "value"
  newNeighbor=(matstate[xNgh+(Ymoins<<DIM)]==value)+(matstate[xNgh+(Yplus<<DIM)]==value)+(matstate[Xmoins+(yNgh<<DIM)]==value)+(matstate[Xplus+(yNgh<<DIM)]==value)+(matstate[Xplus+(Ymoins<<DIM)]==value)+(matstate[Xplus+(Yplus<<DIM)]==value)+(matstate[Xmoins+(Ymoins<<DIM)]==value)+(matstate[Xmoins+(Yplus<<DIM)]==value);

  return newNeighbor;
}


//****************Number of neighbouring new species *******************//
// Following two fxns  determine number of new sp. neighbours and 
// to select one neighbour sp.


//fxn returns number of non-mussel new species in neighbourhood
int XYNeighbour8_newSpecies(int aPos)
{

  xNgh =aPos%SIZE;
  yNgh =(int) (aPos>>DIM); 

  // Define neighborhood coordinates
  Xmoins=xNgh-1;
  Xplus=xNgh+1;
  Ymoins=yNgh-1;
  Yplus=yNgh+1;

  // Apply periodic boundary conditions
  if((xNgh-1)<0) Xmoins=SIZE-1;
  else if((xNgh+1)>(SIZE-1)) Xplus=0;
  if((yNgh-1)<0) Ymoins=SIZE-1;
  else if((yNgh+1)>(SIZE-1)) Yplus=0;

  int count=0;
  int nghbCount[]={matstate[xNgh+(Ymoins<<DIM)], matstate[xNgh+(Yplus<<DIM)],
                   matstate[Xmoins+(yNgh<<DIM)], matstate[Xplus+(yNgh<<DIM)], 
				   matstate[Xplus+(Ymoins<<DIM)], matstate[Xplus+(Yplus<<DIM)], 
				   matstate[Xmoins+(Ymoins<<DIM)], matstate[Xmoins+(Yplus<<DIM)]};

  for(i=0; i<8; i++)
  {
	  if(nghbCount[i]>99)
		  count++;
  }

  
  
  return count;
}



//Function selects and returns one newspecies value from among neighbours
int neighborSpecies(int aPos)
{
  
  xNgh =aPos%SIZE;
  yNgh =(int) (aPos>>DIM); 

  // Define neighborhood coordinates
  Xmoins=xNgh-1;
  Xplus=xNgh+1;
  Ymoins=yNgh-1;
  Yplus=yNgh+1;

  // Apply periodic boundary conditions
  if((xNgh-1)<0) Xmoins=SIZE-1;
  else if((xNgh+1)>(SIZE-1)) Xplus=0;
  if((yNgh-1)<0) Ymoins=SIZE-1;
  else if((yNgh+1)>(SIZE-1)) Yplus=0;

 

  int neighbour=0;
  int count=0;
  int temp_value;
  int nghbCount[]={matstate[xNgh+(Ymoins<<DIM)], matstate[xNgh+(Yplus<<DIM)],
                   matstate[Xmoins+(yNgh<<DIM)], matstate[Xplus+(yNgh<<DIM)],
                   matstate[Xplus+(Ymoins<<DIM)], matstate[Xplus+(Yplus<<DIM)],
		   matstate[Xmoins+(Ymoins<<DIM)], matstate[Xmoins+(Yplus<<DIM)]};


  count= XYNeighbour8_newSpecies(aPos);

  if(count!=0)
  {
	  //int temp[count];
	  vector<int> temp(count);
	  int index=0;
  
	  for (i=0; i<8; i++)
	  {
		  if (nghbCount[i]>99)
                  {
	              temp[index]=nghbCount[i];
		      index++;
                   }
	  }

     temp_value= (int)(randomMT()*count);
     if (temp_value==count)
	     temp_value=count-1;   //check to see actual range of randomMT()

     neighbour= temp[temp_value];
  }
  
  return neighbour;
}


//***************************end************************//

//******function randomly selects/returns a neighbour state if it is a mussel or new sp.*******//
int select_Neighbor(int aPos)
{
  
  xNgh =aPos%SIZE;
  yNgh =(int) (aPos>>DIM); 

  // Define neighborhood coordinates
  Xmoins=xNgh-1;
  Xplus=xNgh+1;
  Ymoins=yNgh-1;
  Yplus=yNgh+1;

  // Apply periodic boundary conditions
  if((xNgh-1)<0) Xmoins=SIZE-1;
  else if((xNgh+1)>(SIZE-1)) Xplus=0;
  if((yNgh-1)<0) Ymoins=SIZE-1;
  else if((yNgh+1)>(SIZE-1)) Yplus=0;

  int neighbour=0;
  
  int nghbArray[]={matstate[xNgh+(Ymoins<<DIM)], matstate[xNgh+(Yplus<<DIM)],
                   matstate[Xmoins+(yNgh<<DIM)], matstate[Xplus+(yNgh<<DIM)],
                   matstate[Xplus+(Ymoins<<DIM)], matstate[Xplus+(Yplus<<DIM)],
		   matstate[Xmoins+(Ymoins<<DIM)], matstate[Xmoins+(Yplus<<DIM)]};


  int index=(int)(8*randomMT());

  if(nghbArray[index]>99 || nghbArray[index]==2)
         neighbour=nghbArray[index];

  return neighbour;



}

//********************end*****************************//


int pos_from_xy(int aX, int anY){
 return aX+(anY<<DIM);
}

int x_from_pos(int aPos){
  return aPos%SIZE; // get x coordinate of "aPos" cell from 1D array
}

int y_from_pos(int aPos){
 return (int) (aPos>>DIM);
 }

int periodic_pos(int aPos){
  int newPos;
  if(aPos<0) return SIZE-aPos;
  else if(aPos>(SIZE-1)) return aPos-SIZE;
  else return aPos;
}






void init_matcover()
{
	for (k=0; k<NSITES;k++)
	{
		matcover[k]= 0;

	}



}







		
//****************************************************************
//The following fxns measure (and rank) non-mussel species abundance,
// ricness and persistence 
//*****************************************************************

//--------updates migrantAbund[] array---------
void migrant_abundance()
{

/*
	int matcopy[NSITES];


	for (k=0; k<NSITES;k++)
	{
		matcopy[k]=matstate[k];

	}

	for (i=0; i<NSITES;i++)
	{
		if (matcopy[i]<100)
		{
			matcopy[i]=0;
		}
	}

*/	
	//initialize migrantAbund[] array to 0
	for(i=0;i<MIGRANTPOOL;i++)
		migrantAbund[i]=0;


	for (i=0;i<NSITES;i++)
	{
		if(matstate[i]>99)
				migrantAbund[matstate[i]-100]++;
	}

	//for(i=0;i<MIGRANTPOOL;i++)
	//{
	//	migrantAbund[i]=speciesCount(i+100);
	//}

}




// ---------measures species richness--------
int richness()
{
        
	int count=0;
	for(i=0;i<MIGRANTPOOL;i++)
	{
		if(migrantAbund[i]>0)
			count++;
	}

	return count;
}



//-----returns number of ind. of a particular sp---------
int speciesCount(int value)
{
	int tempcount=0;
	for(i=0;i<NSITES;i++)
	{
		if(matstate[i]==value)
			tempcount++;
	}

	return tempcount;
}


//****************End of new sp.  measures*******************		    












