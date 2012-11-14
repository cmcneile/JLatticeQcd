//
//  Port Creutz's code to java
//
//

// http://www.cs.geneseo.edu/~baldwin/reference/random.html


//
//


import java.util.Random;


public class generate_gauge {

    public static void main(String[] argv)
    {

	int a , b , c ; 
	int N = Global.GROUP ; 
	int sweeps_between_meas = 5 ;
	int max_sweeps = 5 ; 

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

	for (int iv=0;iv<nlinks;iv++)
	    {
		System.out.println("Matrix " + iv);
		ulinks[iv].printmatrix() ;
	    }

	loop(ulinks,2,2);
	System.exit(0) ;
    

	/* experiment: standard Monte Carlo updating */
	System.out.println("Start Monte Carlo generation of configs");
	for(int iter=0; iter < max_sweeps ; iter++) 
	    {
		int count=0;
		for (int i=0 ; i < sweeps_between_meas ; i++) 
		    {
		    //		    monte(ulinks);
		    count++; 
		}
		renorm(ulinks);
		loop(ulinks,2,2);
	    } 
	




    }


    public static void init()
    {
	int DIM = Global.DIM ; 
	int i , iv ;
	int[] x = new int[DIM] ;
	int N = Global.GROUP ; 

	System.out.println("Initializing lattice");

	nsites=1;
	for(i=0;i<DIM;i++){
	    nsites*=Global.shape[i];
	    //	    if (1 & Global.shape[i]) cleanup(string("bad dimensions"));
	}
	nlinks=DIM*nsites;
	nplaquettes=DIM*(DIM-1)*nsites/2;

	// split lattice into even and odd sites
	vectorlength=nsites/2;

	System.out.println("nlinks = " + nlinks );
	System.out.println("vectorlength = " + vectorlength);

	generator = new Random();

	// 
	//  reserve memory
	//  

	 ulinks=new gaugefield[nlinks];
	 for(i=0 ; i < nlinks ; ++i)
	    {
		ulinks[i] = new gaugefield(N) ;
	    }

	/* set starting links to identity matrix */
	for (iv=0;iv<nlinks;iv++)
	    ulinks[iv].set_unit() ;

	parity=new int[nsites];
	table1=new gaugefield[vectorlength];
	table2=new gaugefield[vectorlength];

	for(i=0 ; i < vectorlength ; ++i)
	    {
		table1[i] = new gaugefield(N);
		table2[i] = new gaugefield(N);
	    }


	mtemp0 =new gaugefield[vectorlength];
	mtemp1 =new gaugefield[vectorlength];
	mtemp2 =new gaugefield[vectorlength];
	mtemp3 =new gaugefield[vectorlength];
	mtemp4 =new gaugefield[vectorlength];

	for(i=0 ; i < vectorlength ; ++i)
	    {
		mtemp0[i] = new gaugefield(N);
		mtemp1[i] = new gaugefield(N);
		mtemp2[i] = new gaugefield(N);
		mtemp3[i] = new gaugefield(N);
		mtemp4[i] = new gaugefield(N);
	    }


	    //	    mtemp[i]=new gaugefield[vectorlength];



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


/**
 * Generate tables of vectorlength random matrices 
 *
 *
 */
    public static void maketable()
    {

	int i,j,iv;
	int GROUP = Global.GROUP ;

	//#define forvector for(iv=0;iv<vectorlength;iv++)
	//#define formatrix for(i=0;i<GROUP;i++)for(j=0;j<GROUP;j++)


	gaugefield temporary1 = new gaugefield(GROUP) ;
	gaugefield temporary2 = new gaugefield(GROUP) ;
	
	for(iv=0;iv<vectorlength;iv++)
	    {
		/* bias towards the identity */
		temporary1.set_constant(Global.beta/GROUP, 0.0);
		temporary2.set_constant(Global.beta/GROUP, 0.0) ;
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
	//	ranmat(mtemp[0]);
	mtemp0 = ranmat();

	/* multiply table 2 by a into table 1 for trial change */
	//	table1 =  vprod(table2 , mtemp[0]);
	table1 =  vprod(table2 , mtemp0 );

	/* metropolis select new table 2 */
	sold = vtrace(table2 );
	snew =  vtrace(table1);
	metro(table2,table1,6*Global.beta/Global.GROUP);  

	/* switch table 1 and 2 */
	table1 =  vcopy(table2);
	table2 = vcopy(mtemp0); 

	vgroup(table1);

	return;
    }


  /* 
     set g3 to the matrix product of g1 and g2, vectorlength times 
  */    

    public static gaugefield[] vcopy(gaugefield[] g1) 
    {
	gaugefield[] g3 ; 
	g3 = new gaugefield[vectorlength] ;
	
	for(int iv=0;iv<vectorlength;iv++)
	    g3[iv]=g1[iv].copy() ;
	
	return g3 ;
    }

  /* matrix sum of g1 and g2 to g3, vectorlength times */
    public static gaugefield[] vsum(gaugefield[] g1,gaugefield[] g2) 
    {
	gaugefield[] g3 ;
	g3 = new gaugefield[vectorlength] ;
	

	int i,j,iv;
	/*   slightly faster writing this out over:
	     forvector
	     g3[iv]=g1[iv]+g2[iv];
	*/  
	for(i=0;i<Global.GROUP;i++)
	    for(j=0;j<Global.GROUP;j++)
		for(iv=0;iv<vectorlength;iv++) {
		    g3[iv].real[i][j]=g1[iv].real[i][j]+g2[iv].real[i][j];
		    g3[iv].imag[i][j]=g1[iv].imag[i][j]+g2[iv].imag[i][j];
		}

	return g3 ; 
    }


  /* 
     randomly shift table1, randomly invert, and put in g 
  */
  public static gaugefield[] ranmat() 
    {
	gaugefield[] g ; 
	g = new gaugefield[vectorlength] ;

	// index=(int) (vectorlength*drand48());
	int index = generator.nextInt(vectorlength) + 1 ;

	for(int iv=0;iv<vectorlength;iv++)
	    {
		if (index>=vectorlength) index-=vectorlength;
		
		double rr = generator.nextDouble();
		
		if( rr < 0.5 )
		    g[iv]=table1[index].copy() ;
		else
		    g[iv]=table1[index].conjugate();

		index++;
	    }
	
	return g;
    }



  public static double[] vtrace(gaugefield[] g1) {
    double[] g3 ; 
    g3 = new double[vectorlength] ;

    for(int iv=0;iv<vectorlength;iv++)
	g3[iv]=g1[iv].trace_re() ;


  return g3 ;
}



  /* 
     set g3 to the matrix product of g1 and g2, vectorlength times 
  */    

  public static  gaugefield[] vprod(gaugefield[] g1 , gaugefield[] g2) 
{
    gaugefield[] g3 ; 
    g3 = new gaugefield[vectorlength] ;

    for(int iv=0;iv<vectorlength;iv++)
	g3[iv]=g1[iv].prod(g2[iv]) ;


  return g3 ;
}


    /* the basic metropolis update */
    public static double monte(gaugefield[]  lattice) 
    {
	int DIM = Global.DIM ;
	int HITS = Global.HITS ;
	int GROUP = Global.GROUP ;
	
	double stot,acc,eds;
	int iv,iacc,color,link,hit;

	/* update table */
	vtable(); 
	stot=eds=0.0;
	iacc=0;
	/* loop over checkerboard colors */

	for (color=0;color<2;color++) {
	    /* loop over link dirs */
	    for (link=0;link<DIM;link++) {
		/* get neighborhood */
		mtemp4 = staple(lattice,color,link); 
		/* get old link and calculate action */
		mtemp0 = getlinks(lattice,color,link);
		sold = vtprod(mtemp0,mtemp4);
		/* loop over hits */
		for (hit=0;hit<HITS;hit++) {
		    /* get random matrices */
		    mtemp1 = ranmat() ;

		    /* find trial element and new action */
		    mtemp2 = vprod(mtemp0,mtemp1);
		    snew = vtprod(mtemp2,mtemp4);
		    eds += metro(mtemp0,mtemp2,
				 Global.beta/(1.*Global.GROUP)); 

		    /* metropolis step */
		    for(iv=0;iv<vectorlength;iv++) {
			iacc=iacc+accepted[iv];
			stot=stot+sold[iv];
		    }  
		}
		lattice = savelinks(mtemp0,color,link); 
	    }
	}
	stot=stot/(.5*DIM*(DIM-1)*nlinks*GROUP*HITS);
	acc=iacc/(1.*nlinks*HITS);
	eds=eds/(2.*DIM*HITS);


	/* eds should fluctuate about unity when in equilibrium */
	System.out.printf("stot=%f, acc=%f, eds=%f\n",stot,acc,eds);


	return stot;
}


  /* This subroutine calculates a vector of matrices interacting with
     with links using Wilson action.  The lattice is in lat and the
     result is placed in st.  The first three matrix vectors mtemp[0],
     mtemp[1], and mtemp[2], are used; so st should not be there and these
     shouldn't be used until after staple is done.  
     links and sites labeled as

     2--link2--x
     link3     link1
     0--link --1
     link6     link4
     5--link5--4
  */


public static gaugefield[] staple(gaugefield[] lat,int site,int link) 
{
  int iv,link1,site1,site2,site4,site5;
  gaugefield[] st ; // MoreWork need to reserve memory

  st = new gaugefield[vectorlength] ; 

  for(iv=0;iv<vectorlength;iv++)
      st[iv] = new gaugefield(Global.GROUP)  ;

  site1=ishift(site,link,1);
  /* loop over planes */
  for (link1=0;link1<Global.DIM;link1++)
    if (link1!=link) {
      site2=ishift(site ,link1, 1);
      site4=ishift(site1,link1,-1);
      site5=ishift(site ,link1,-1);
      /* top of staple */
      mtemp0 = getlinks(lat,site1,link1);
      mtemp1 = getconjugate(lat,site2,link);
      mtemp2 = vprod(mtemp0,mtemp1);
      mtemp0 = getconjugate(lat,site,link1);
      mtemp1 = vprod(mtemp2,mtemp0);
      st = vsum(st,mtemp1);
      /* bottom of staple */
      mtemp0 = getconjugate(lat,site4,link1);
      mtemp1 = getconjugate(lat,site5,link );
      mtemp2 = vprod(mtemp0,mtemp1);
      mtemp0 = getlinks(lat,site5,link1);
      mtemp1 = vprod(mtemp2,mtemp0);
      st = vsum(st,mtemp1);
    }


  return st ;
}


    /* gather conjugate links into vector g */
    public static gaugefield[] getconjugate(gaugefield[] lat,int site,int link)
    {
	gaugefield[] g = new gaugefield[vectorlength] ; 
	int iv,shift;
	makeindex(site,myindex);
	shift=nsites*link;

	for(iv=0;iv<vectorlength;iv++)
	    g[iv]=lat[myindex[iv]+shift].conjugate();


  return g ;
}

    /* generates a set of site labels starting at n for gathering links */
    /* loop over even parity sites and gather with shift n from them */


public static void makeindex(int n,int[] ind){

  int iv,site;
  int[] x = new int[Global.DIM] ;

  split(x,n);
  site=0;

  for(iv=0;iv<vectorlength;iv++)
      {
    while (parity[site] != 0 ) site++;
    ind[iv]=vshift(site,x);
    site++;
  }
  return;
}



    /**
     * Shift a site n by vector x
     *
     * Longer description. If there were any, it would be    [2]
     * here.
     *
     * @param  n  starting site
     * @param  x[] shift direction
     * @return site (integer)
     */


    public static int vshift(int n, int[] x)
    {
	int i ;
	int[] y = new int[Global.DIM];
	
	split(y,n);
	for(i=0;i<Global.DIM;i++){
	    if (x[i] != 0){
		y[i]+=x[i];
		while (y[i]>=Global.shape[i])
		    y[i]-=Global.shape[i];
		while (y[i]<0)
		    y[i]+= Global.shape[i];
	    }
	}
	return siteindex(y);

}


  /* utilities for implementing periodic boundaries */

  /* gives a unique index to site located at x[DIM] */
    public static int siteindex(int[] x){

  int i,result=0;

  for (i=0;i<Global.DIM;i++)
    result += shift[i]*x[i];

  return result;
}



    /**
     * returns index of a site shifted dist in direction dir
     *
     * Longer description. If there were any, it would be    [2]
     * here.
     *
     * @param  n
     * @param  direction   -- direction that the sift is in
     * @param  dist
     * @return site (integer)
     */

    public static int ishift( int n, int dir, int dist)
    {
	
	int i ;
	int[] x = new int[Global.DIM];
	for(i=0;i<Global.DIM;i++)
	    x[i]=0;

	x[dir]=dist;
	return vshift(n,x);
}

  /* gather same color links into vector g starting at site */
public static gaugefield[]  getlinks(gaugefield[] lattice,int site,int link)
{
    gaugefield[] g = new gaugefield[vectorlength] ;
    /** more work memory  ***/

    int iv,shift;
    makeindex(site,myindex);
    shift=nsites*link;

    for(iv=0;iv<vectorlength;iv++)
	g[iv]=lattice[myindex[iv]+shift];


  return g ;
}

/* real trace of product g1 and g2 to s, vectorlength times */

    public static double[] vtprod(gaugefield[] g1, gaugefield[] g2)
{
    double[] s = new double[vectorlength] ;

    for(int iv=0;iv<vectorlength;iv++)
	s[iv]=0.0;


    for(int iv=0;iv<vectorlength;iv++)
	for(int i=0;i<Global.GROUP;i++)
	    for(int j=0;j<Global.GROUP;j++)
		s[iv]+=g1[iv].real[i][j]*g2[iv].real[j][i]
		    -g1[iv].imag[i][j]*g2[iv].imag[j][i];   
    
    return s ;
}
    

  /* accept new for old using metropolis algorithm
     return average exponential of action change, this should fluctuate
     about unity when in equilibrium
     bias multiplies actions in exponential (i.e. beta/GROUP)  
     accepted changes returned in accepted
     actions passed in global variables sold and snew */


    public static double metro(gaugefield[] old,gaugefield[] trial,double bias) 
    {
	int iv;
	double expdeltas=0.0,temp;
	

	for(iv=0;iv<vectorlength;iv++)
	    {
		temp=  Math.exp((bias*(snew[iv]-sold[iv])));
		expdeltas=expdeltas+temp;
		
		// this is a big hack from the c++ cpde
		if(generator.nextDouble()  < temp)
		    accepted[iv] = 1 ;
		else
		    accepted[iv] = 0 ;
		
	    }
	/* Accept changes */

	for(iv=0;iv<vectorlength;iv++)
	    if (accepted[iv] != 0 ) 
		{
		    sold[iv]=snew[iv];
		    old[iv]=trial[iv];
		}


  return expdeltas/vectorlength;
}   


  /* scatter alternate links from vector g */

public static gaugefield[] savelinks(gaugefield[]  g,int site,int link)
{
    gaugefield[] lattice = new gaugefield[vectorlength] ;

    makeindex(site,myindex);
    int shift=nsites*link;

    for(int iv=0;iv<vectorlength;iv++)
      lattice[myindex[iv]+shift] = g[iv].copy();

    return lattice ;
}


/**
 * calculate rectangular wilson loops 
 *
 * Compute Wilson loop with size x * y
 *
 * @param gaugefield -- gauge configuration
 * @param x -- one dimension of Wilson loop
 * @param y -- the other dimension of Wilson loop
 */
    public static double loop(gaugefield[] u, int x, int y)
    {
	int i,color,link1,link2,iv,corner,count=0;
	int DIM = Global.DIM ; 

	double result=0.;
	for(color=0;color<2;color++)
	    for (link1=0;link1<DIM;link1++)
		{

		    //		for (link2=(x==y)*(link1+1);link2<DIM;link2++)
		    int link2_start ; 
		    if( x==y )
			link2_start = link1+1 ; 
		    else
			link2_start = 0 ; 

		for (link2=link2_start ;link2<DIM;link2++)
		    if (link1 != link2){
			count++;
			corner=ishift(color,link1,x);
			corner=ishift(corner,link2,y);
			for(iv=0;iv<vectorlength;iv++)
			    {
				mtemp0[iv].set_unit() ;
				mtemp1[iv].set_unit() ;
				mtemp2[iv].set_unit() ;
				mtemp3[iv].set_unit() ;
			    }
			for(i=0;i<x;i++)
			    {
				mtemp4 =  getlinks(u,ishift(color,link1,i),link1);
				mtemp0 = vprod(mtemp0,mtemp4);
				mtemp4 = getconjugate(u,ishift(corner,link1,-i-1),link1);
				mtemp2 = vprod(mtemp2,mtemp4);
			    }

			for(i=0;i<y;i++)
			    {
				mtemp4 = getlinks(u,ishift(corner,link2,i-y),link2);
				mtemp1 = vprod(mtemp1,mtemp4);
				mtemp4 = getconjugate(u,ishift(color,link2,y-i-1),link2);
				mtemp3 = vprod(mtemp3,mtemp4);
			    }
			mtemp0 = vprod(mtemp0,mtemp1);
			mtemp0 = vprod(mtemp0,mtemp2);
			sold = vtprod(mtemp0,mtemp3);
	  
			for(iv=0;iv<vectorlength;iv++)
			    result +=sold[iv];

		    } // link1 != link2
		}

	result=result/(Global.GROUP*vectorlength*count);

	System.out.printf("W[%d,%d] = %g\n",x,y,result);
	return result;
    }

/**
 * Reunitarise gauge configuratiom.
 *
 * project whole lattice into group; call from time to time 
 *    to keep down drift from floating point errors 
 *
 *
 * @param  gaugefield -- gauge configuration
 *
 */


public static void renorm(gaugefield[]  l) 
{ 
  int iv,octant,link;
  /* loop over lattice octants */
  for (octant=0;octant<2*Global.DIM;octant++) 
      {
	  link=octant*vectorlength;
	  for(iv=0;iv<vectorlength;iv++)
	      l[link+iv].project();
      }
  return;
}


    //
    // variables in the class
    //

    public static int nsites ;
    public static int nlinks ; 
    public static int nplaquettes ;
    public static int vectorlength; 


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

    public static gaugefield[] mtemp0 ;
    public static gaugefield[] mtemp1 ;
    public static gaugefield[] mtemp2 ;
    public static gaugefield[] mtemp3 ;
    public static gaugefield[] mtemp4 ;

    // Random number generator 
    private static Random generator ;
}
