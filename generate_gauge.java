//
//  Port Creutz's code to java
//
//

// http://www.cs.geneseo.edu/~baldwin/reference/random.html

import java.util.Random;


public class generate_gauge {

    public static void main(String[] argv)
    {

	int a , b , c ; 
	int N = Global.GROUP ; 
	int[] shape = {8, 8 ,8 , 8 } ; 

	beta = 6.0 ; 
	
	System.out.println("Pure gauge simulation SU(" + N + ")");
	init() ;

	System.out.printf("lattice size %d " ,  Global.shape[0] );
	for (int i=1;i< Global.DIM;i++)
	    System.out.printf(" by %d ", Global.shape[i]);

	System.out.printf("\n");

	//    System.out.println("\n vectorlength = " + vectorlength);
	System.out.println("group=SU(" + N +   " )   beta = " + 
			   Global.beta);
	System.out.println("-----------------");
    

	//	gaugefield  u(N) ; 

    }


    public static void init()
    {
	int DIM = Global.DIM ; 
	int i , iv ;
	int[] x = new int[DIM] ;
	int N = Global.GROUP ; 

	System.out.println("Initializing lattice");

	for(i=0;i<DIM;i++){
	    nsites*=Global.shape[i];
	    //    if (1 & Global.shape[i]) cleanup(string("bad dimensions"));
	}
	nlinks=DIM*nsites;
	nplaquettes=DIM*(DIM-1)*nsites/2;
	vectorlength=nsites/2;

	Random generator = new Random();


	// 
	//  reserve memory
	//  

	 ulinks=new gaugefield[nlinks];
	for(i=0 ; i < nlinks ; ++i)
	    {
		ulinks[i] = new gaugefield(N) ;
	    }

	parity=new int[nsites];
	table1=new gaugefield[vectorlength];
	table2=new gaugefield[vectorlength];
	for(i=0 ; i < vectorlength ; ++i)
	    {
		table1=new gaugefield(N);
		table2=new gaugefield(N);
	    }

	// for (i=0;i<5;i++)
	// mtemp[i]=new matrix[vectorlength];


	sold=new double[vectorlength];
	snew=new double[vectorlength];
	accepted=new int[vectorlength];
	myindex=new int[vectorlength];
	myindex2=new int[vectorlength];

	/* initialize shift array for locating links */
	shift =new int[DIM];
	shift[0]=1;
	for (i=1;i<DIM;i++)
	    shift[i]=shift[i-1]*Global.shape[i-1];

	/* set starting links to identity matrix */
	//  for (iv=0;iv<nlinks;iv++)
	//  ulinks[iv]=1.0;
	
	/* set parity matrix for sites */
	for (iv=0;iv<nsites;iv++)
	    {
		split(x,iv);
		parity[iv]=0;
		for (i=0;i<DIM;i++)
		parity[iv] ^= x[i];
		parity[iv] &= 1;
	    }

  maketable();
  System.out.println("initialization done\n");
  

    }



    public static void maketable(){
  /* generate tables of vectorlength random matrices */
  int i,j,iv;
  int GROUP = Global.GROUP ;

  //#define forvector for(iv=0;iv<vectorlength;iv++)
  //#define formatrix for(i=0;i<GROUP;i++)for(j=0;j<GROUP;j++)

  Random generator = new Random();

  gaugefield temporary1 = new gaugefield(GROUP) ;
  gaugefield temporary2 = new gaugefield(GROUP) ;

  for(iv=0;iv<vectorlength;iv++)
      {
	  /* bias towards the identity */
	  temporary1.set_constant(beta/GROUP, 0.0);
	  temporary2.set_constant(beta/GROUP, 0.0) ;
	      for(i=0;i<GROUP;i++)
		  for(j=0;j<GROUP;j++)
		      {
			  temporary1.real[i][j] += 
			      generator.nextDouble() -0.5;
			  temporary1.imag[i][j] += 
			      generator.nextDouble() -0.5;
			  temporary2.real[i][j] += 
			      generator.nextDouble() -0.5;
			  temporary2.imag[i][j] += 
			      generator.nextDouble() -0.5;
		      } 
	  table1[iv]=temporary1;
	  table2[iv]=temporary2;
      }
  /* make into group elements */
  vgroup(table1);
  vgroup(table2);
  /* update table a few times */
  for (i=0;i<50;i++)
    vtable();
  return;
}

    /* 
       subroutine to make group elements out of 
       vectorlength matrices g 
    */
    public static void vgroup(gaugefield[] g)
    {
	int iv;
	for(iv=0;iv<vectorlength;iv++)
	    g[iv].project();
	return;
    }





    //
    //  Basic end of program
    //
    public static void cleanup(String msg )
    {
	System.out.println("Error: " + msg);
	System.exit(1) ;
    }


    /* splits a site index into coordinates 
       assume s in valid range 0<=s<nsites 
       I (Creutz) think this is faster than using mods, 
       but this should be tested 
    */
 

    public static void split(int[] x, int s)
    {
	int DIM = Global.DIM ; 

	int i;
	if (s<0 || s>=nsites) cleanup("bad split");
	for(i=DIM-1;i>0;i--){
	    x[i]=0;
	    while (s>=shift[i]){
		s-=shift[i];
		x[i]++;
	    }
	}
	x[0]=s;
	return;
    }


    /* update matrix table */
    public static void vtable() 
    {/* shuffle table 1  into a */
	/* the random inversion from ranmat is important! */
	//  ranmat(mtemp[0]);
	/* multiply table 2 by a into table 1 for trial change */
	// vprod(table2,mtemp[0],table1);
	/* metropolis select new table 2 */
	// vtrace(table2,sold);
	// vtrace(table1,snew);
	// metro(table2,table1,6*beta/GROUP);  
	/* switch table 1 and 2 */
	// vcopy(table2,table1);
	// vcopy(mtemp[0],table2);
	// vgroup(table1);
	return;
    }
    


    //
    // variables in the class
    //

    public static int nsites ;
    public static int nlinks ; 
    public static int nplaquettes ;
    public static int vectorlength; 

    public double beta ;

    public static int[] accepted ;
    public static int[] myindex ;
    public static int[] myindex2 ;
    public static int[] parity ;

    public static double[] sold ;
    public static double[] snew ;

    public static int[] shift ;

    public static gaugefield[] ulinks; /* for the main lattice */

    public static gaugefield[] table1 ;
    public static gaugefield[] table2 ;

}
