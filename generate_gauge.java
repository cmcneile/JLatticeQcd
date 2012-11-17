/**
 * @author      Craig McNeile
 * @version     0.5      
 * @since       2012-11-1
 *
 *  This is a rewrite of a C++ code written by Mike Creutz.
 *
 */




import java.util.Random;

import java.io.IOException;
import org.w3c.dom.*;
import org.xml.sax.SAXException;
import javax.xml.parsers.*;
import javax.xml.xpath.*;


public class generate_gauge {

    public static void main(String[] argv)
	throws ParserConfigurationException, SAXException, 
	       IOException, XPathExpressionException
    {

	//	int a , b , c ; 
	int sweeps_between_meas = 5 ;
	int max_sweeps = 50 ; 


	readparam param_file = new readparam()  ;

	String xpath_want = "//pure_gauge_params/algorithm_param/max_sweeps/text()" ;
	max_sweeps = param_file.read_int(xpath_want) ;

	N = Global.GROUP ; 

	System.out.println("Pure gauge simulation of SU(" + N + ") theory");
	init() ;

	System.out.printf("Lattice size %d" ,  Global.shape[0] );
	for (int i=1;i< Global.DIM;i++)
	    System.out.printf("x%d", Global.shape[i]);

	System.out.printf("\n");

	System.out.println("Gauge group=SU(" + N +   ")   beta = " + 
			   Global.beta);
	System.out.println("-----------------");

	//	dump_matrix(ulinks); 
	// System.exit(0) ;

	loop(ulinks,1,1);
	

    

	/* Standard Monte Carlo updating */
	System.out.println("Start Monte Carlo generation of configs");
	for(int iter=0; iter < max_sweeps ; iter++) 
	    {
		int count=0;
		for (int i=0 ; i < sweeps_between_meas ; i++) 
		    {
			//			System.out.println("DEBUG Sweep number " + i);
			monte(ulinks);

			//			check_unitarity(ulinks) ;
			// System.exit(0) ;

			count++; 
		    }
		check_unitarity_norm(ulinks) ; // System.exit(0) ; 
		renorm(ulinks);
		loop(ulinks,1,1);
		check_unitarity_norm(ulinks) ; // System.exit(0) ; 

	    } 
	




    }




/**
 * Reserve memory and initialize mapping arrays.
 *
 */
    public static void init()
    {
	int DIM = Global.DIM ; 
	int i , iv ;
	int[] x = new int[DIM] ;

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


	//
	// arrays with nlinks members
	//

	 ulinks=new gaugefield[nlinks];
	 for(i=0 ; i < nlinks ; ++i)
	    {
		ulinks[i] = new gaugefield(N) ;
	    }

	 /* set starting links to identity matrix */
	 for (iv=0;iv<nlinks;iv++)
	     ulinks[iv].set_unit() ;

	 //
	 // arrays with nsites members
	 //

	parity=new int[nsites];

	//
	//  arrays with vectorlength meebers
	//

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


	sold=new double[vectorlength];
	snew=new double[vectorlength];

	accepted=new int[vectorlength];

	myindex=new int[vectorlength];

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

	// debug_init(shift) ; System.exit(0) ;

	maketable();
	System.out.println("Initialisation done\n");
  

    }


/**
 * Generate tables of vectorlength random matrices 
 *
 * The table1 and table2 are created.
 *
 */
    public static void maketable()
    {
	gaugefield temporary1 = new gaugefield(N) ;
	gaugefield temporary2 = new gaugefield(N) ;
	
	for(int iv=0;iv<vectorlength;iv++)
	    {
		/* bias towards the identity */
		temporary1.set_diag(Global.beta/N, 0.0);
		temporary2.set_diag(Global.beta/N, 0.0) ;
		for(int i=0; i < N ; i++)
		    for(int j=0; j < N ; j++)
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
		table1[iv]=temporary1.copy() ;
		table2[iv]=temporary2.copy() ;
	    }

	/* make into group elements */
	vgroup(table1);
	vgroup(table2);

	check_unitarity_norm(table2) ; // System.exit(0) ;

	//System.out.println("DEBUG Table 1");
	//dump_matrix(table1) ;

	//System.out.println("DEBUG Table 2");
	//dump_matrix(table2) ;

	//	System.exit(0) ;
	
	/* update table a few times */
	for (int i=0;i<50;i++)
	    {
		//		System.out.println("Updating table " + i);
		vtable();
	    }
	return;
    }




/**
 * Basic end of program 
 *
 *
 * @param  string msg 
 */
    public static void cleanup(String msg )
    {
	System.out.println("Error: " + msg);
	System.exit(1) ;
    }




/**
 *
 *    splits a site index into coordinates 
 *      assume s in valid range 0<=s<nsites 
 *      I (Creutz) think this is faster than using mods, 
 *      but this should be tested 
 *
 *
 * @param  int s
 * @param  int[] 
 */
    public static void split(int[] x, int s)
    {
	int DIM = Global.DIM ; 

	if (s<0 || s>=nsites) cleanup("bad split");

	for(int i=DIM-1 ; i>0 ;i--)
	    {
		x[i]=0;
		while (s>=shift[i])
		    {
			s-=shift[i];
			x[i]++;
		    }
	    }
	x[0]=s;
	return;
    }


   


/**
 * Randomly permutate the table1 and table2 tables
 *
 * The table1 and table2 are modified.
 *  update matrix table
 *
 */
    public static void vtable() 
    {
	mtemp0 = ranmat();

	/* multiply table 2 by a into table 1 for trial change */
	table1 =  vprod(table2 , mtemp0 );

	//System.out.println("DEBUG Table 1"); 
	//dump_matrix(table1) ; 

	//System.out.println("DEBUG Table 2"); 
	//dump_matrix(table2) ; 
	//System.exit(0) ;  

	/* metropolis select new table 2 */
	sold =  vtrace(table2);
	snew =  vtrace(table1);

	// for(int ii=0 ; ii < snew.length ; ++ii)
	//    {
	//	System.out.println("DEBUG old,new " + ii + " " + sold[ii] + " " + snew[ii] ); 
	//    } 

	//	System.exit(0) ;  
	metro(table2,table1,6*Global.beta/N);  

	/* switch table 1 and 2 */
	table1 = vcopy(table2);
	table2 = vcopy(mtemp0); 

	vgroup(table1);

	return;
    }



/**
 * Randomly shift table1, randomly invert, and put in g 
 *
 * @return gaugefield[vectorlength] g1 -- randomly permuted table1
 */


  public static gaugefield[] ranmat() 
    {
	gaugefield[] g ; 
	g = new gaugefield[vectorlength] ;
	for(int i=0 ; i < vectorlength ; ++i)
	    {
		g[i] = new gaugefield(N);
	    }	


	//	int index = generator.nextInt(vectorlength) + 1 ;
	int index = generator.nextInt(vectorlength)  ;

	//		System.out.println("DEBUG random index = " + index);

	/* the random inversion from ranmat is important! */

	for(int iv=0;iv<vectorlength;iv++)
	    {
		if (index>=vectorlength) index-=vectorlength;
		
		double rr = generator.nextDouble();
		//		System.out.println("DEBUG random = " + rr);
		
		if( rr < 0.5 )
		    g[iv]=table1[index].copy() ;
		else
		    g[iv]=table1[index].conjugate();

		index++;
	    }
	
	return g;
    }

    /*****************************************

      Start of vector routines 

     *****************************************/


/**
 * Project an array of matrics to the SU(3) group.
 *
 *
 *
 * @param  gaugefield[vectorlength] g -- vector of gaugefields
 */


    public static void vgroup(gaugefield[] g)
    {
	for(int iv=0;iv<vectorlength;iv++)
	    g[iv].project();

	return;
    }


/**
 * At each vector position compute the real trace of gaugefield
 *
 *
 * @param  gaugefield[vectorlength] g1
 * @return double[vectorlength]          -- real trace of gaugefields
 */

    public static double[] vtrace(gaugefield[] g1) 
    {
	double[] spur ; 
	spur = new double[vectorlength] ;

	for(int iv=0;iv<vectorlength;iv++)
	    spur[iv]=g1[iv].trace_re() ;
    

	return spur ;
    }


/**
 * Vector copy of gauge matrices
 *
 *
 * @param  gaugefield[vectorlength] ginput
 * @return gaugefield[vectorlength] gcopy  -- copy of ginput 
 */
    public static gaugefield[] vcopy(gaugefield[] ginput) 
    {
	gaugefield[] gcopy ; 
	gcopy = new gaugefield[vectorlength] ;
	for(int i=0 ; i < vectorlength ; ++i)
	    {
		gcopy[i] = new gaugefield(N);
	    }	

	for(int iv=0;iv<vectorlength;iv++)
	    gcopy[iv]=ginput[iv].copy() ;
	
	return gcopy ;
    }





/**
 * Matrix sum for a vector of matrices.
 *
 * Longer description. If there were any, it would be    [2]
 * here.
 *
 *
 * @param  gaugefield[vectorlength] g1
 * @param  gaugefield[vectorlength] g2
 * @return double[vectorlength]          -- matrix product g1 + g2 
 */
    public static gaugefield[] vsum(gaugefield[] g1,gaugefield[] g2) 
    {
	gaugefield[] g3 ;
	g3 = new gaugefield[vectorlength] ;
	for(int i=0 ; i < vectorlength ; ++i)
	    {
		g3[i] = new gaugefield(N);
	    }	

	for(int i=0;i< N;i++)
	    for(int j=0;j< N;j++)
		for(int iv=0;iv<vectorlength;iv++) 
		    {
			g3[iv].real[i][j]=g1[iv].real[i][j]+g2[iv].real[i][j];
			g3[iv].imag[i][j]=g1[iv].imag[i][j]+g2[iv].imag[i][j];
		    }

	return g3 ; 
    }




/**
 * Write to standard output the matrices of gauge field
 *
 * Longer description. If there were any, it would be    [2]
 * here.
 *

 * @param  gaugefield[vectorlength] g1
 * @return double[vectorlength]          -- trace of gaugefields
 */

    public static void dump_matrix(gaugefield[] g) 
    {
	for(int iv=0; iv < g.length ;iv++)
	    {
		System.out.println("Matrix " + iv);
		g[iv].printmatrix() ;
	    }
    }




/**
 *    set g3 to the matrix product of g1 and g2, vectorlength times 
 *
 * Longer description. If there were any, it would be    [2]
 * here.
 *
 * @param  gaugefield[vectorlength] g2
 * @param  gaugefield[vectorlength] g1
 * @return double[vectorlength]          -- g3 = g1 * g2
 */
    public static  gaugefield[] vprod(gaugefield[] g1 , gaugefield[] g2) 
    {
	gaugefield[] g3 ; 
	g3 = new gaugefield[vectorlength] ;
	for(int i=0 ; i < vectorlength ; ++i)
	    {
		g3[i] = new gaugefield(N);
	    }	
	
	for(int iv=0;iv<vectorlength;iv++)
	    g3[iv]=g1[iv].prod(g2[iv]) ;

	return g3 ;
    }



/**
 * gather conjugate links into vector g 
 *
 * Longer description. If there were any, it would be    [2]
 * here.
 *
 * @param  int link -- direction
 * @param  int site
 * @param  gaugefield[] lat
 * @return gaugefield[vectorlength]
 */
    public static gaugefield[] getconjugate(gaugefield[] lat,int site,int link)
    {
	gaugefield[] g = new gaugefield[vectorlength] ; 
	for(int i=0 ; i < vectorlength ; ++i)
	    {
		g[i] = new gaugefield(N);
	    }	


	makeindex(site,myindex);
	int shift=nsites*link;

	for(int iv=0;iv<vectorlength;iv++)
	    g[iv]=lat[myindex[iv]+shift].conjugate();


	return g ;
    }



/**
 * real trace of product g1 and g2 to s, vectorlength times 
 *
 * Longer description. If there were any, it would be    [2]
 * here.
 *
 * @param  gaugefield[] g2 --
 * @param  gaugefield[] g1 -- 
 * @return double []       --
 */

    public static double[] vtprod(gaugefield[] g1, gaugefield[] g2)
    {
	double[] s = new double[vectorlength] ;
	
	for(int iv=0;iv<vectorlength;iv++)
	    s[iv]=0.0;


	for(int iv=0;iv<vectorlength;iv++)
	    for(int i=0;i<N;i++)
		for(int j=0;j<N;j++)
		    s[iv]+=g1[iv].real[i][j]*g2[iv].real[j][i]
			-g1[iv].imag[i][j]*g2[iv].imag[j][i];   
	
	return s ;
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
	/* loop over lattice octants. The 2 is for 
	   two parities.
	 */
	for (int octant=0;octant<2*Global.DIM;octant++) 
	    {
		int link=octant*vectorlength;
		for(int iv=0;iv<vectorlength;iv++)
		    l[link+iv].project();
	    }
	return;
    }





    /**
     * Check that the gauge configuration is unitarity
     *
     *
     *
     * @param  gaugefield[nlinks] -- gauge configuration
     *
     */
    

    public static void check_unitarity(gaugefield[]  l) 
    { 
	/* loop over lattice octants. The 2 is for 
	   two parities.
	 */
	for (int octant=0;octant<2*Global.DIM;octant++) 
	    {
		int link=octant*vectorlength;
		for(int iv=0;iv<vectorlength;iv++)
		    {
			System.out.printf("Position[%d,%d]\n" ,link, iv   );
			l[link+iv].check_unitarity();
		    }

	    }
	return;
    }





    /**
     * Check that the gauge configuration is unitarity
     *
     * This works for a vector of any length.
     *
     * @param  gaugefield[] -- gauge configuration
     *
     */
    

    public static void check_unitarity_norm(gaugefield[]  l) 
    { 
	double norm_max = 0.0 ; 
	double dim = l.length ;

	
	for(int iv=0;iv < dim ; iv++)
	    {
		double norm = l[iv].check_unitarity_norm();
		if( norm > norm_max ) norm_max = norm ; 
	    }
	System.out.printf("Maximum  norm = %g\n" , norm_max); 

	return;
    }







    /******************************************

       MonteCarlo routines.

     ****************************************/


/**
 * The basic metropolis update
 *
 * Longer description. If there were any, it would be    [2]
 * here.
 *
 * @param  gaugefield[] lattice configurations (direction and spacetime)
 */
    public static double monte(gaugefield[]  lattice) 
    {
	int DIM = Global.DIM ;
	int HITS = Global.HITS ;
	
	double stot,acc,eds;

	/* update table */
	vtable(); 

	stot=eds=0.0;
	int iacc = 0;

	/* loop over checkerboard colors */

	for (int color=0;color<2;color++) 
	    {
		/* loop over link dirs */
		for (int link=0;link<DIM;link++) 
		    {
			/* get neighborhood */
			mtemp4 = staple(lattice,color,link); 

			/* get old link and calculate action */
			mtemp0 = getlinks(lattice,color,link);
			sold = vtprod(mtemp0,mtemp4);

			/* loop over hits */
			for (int hit=0;hit<HITS;hit++) 
			    {
				/* get random matrices */
				mtemp1 = ranmat() ;

				/* find trial element and new action */
				mtemp2 = vprod(mtemp0,mtemp1);
				snew = vtprod(mtemp2,mtemp4);
				eds += metro(mtemp0,mtemp2,
					     Global.beta/(1.*N)); 

				/* metropolis step */
				for(int iv=0;iv<vectorlength;iv++) 
				    {
					iacc=iacc+accepted[iv];
					stot=stot+sold[iv];
				    }  
			    }
			savelinks(lattice, mtemp0,color,link); 
		    }
	    }

	stot=stot/(.5*DIM*(DIM-1)*nlinks*N*HITS);
	acc=iacc/(1.*nlinks*HITS);
	eds=eds/(2.*DIM*HITS);


	/* eds should fluctuate about unity when in equilibrium */
	System.out.printf("stot=%f, acc=%f, eds=%f\n",stot,acc,eds);


	return stot;
}



/**
 * Compute the staple of gauge fields
 *
 *
 * This subroutine calculates a vector of matrices interacting with
 *   with links using Wilson action.  The lattice is in lat and the
 *    result is placed in st.  The first three matrix vectors mtemp[0],
 *    mtemp[1], and mtemp[2], are used; so st should not be there and these
 *    shouldn't be used until after staple is done.  
 *    links and sites labeled as
 *
 *    2--link2--x
 *    link3     link1
 *    0--link --1
 *    link6     link4
 *     5--link5--4
 *
 *
 * @param  int link -- 
 * @param  int site -- 
 * @param  gaugefield[] latt configurations (direction and spacetime)
 * @return gaugefield[] 
 */



public static gaugefield[] staple(gaugefield[] lat,int site,int link) 
{
  int site1,site2,site4,site5;
  gaugefield[] st ; 

  st = new gaugefield[vectorlength] ; 
  for(int i=0 ; i < vectorlength ; ++i)
      {
	  st[i] = new gaugefield(N);
      }	

  site1=ishift(site,link,1);

  /* loop over planes */
  for (int link1=0;link1<Global.DIM;link1++)
    if (link1!=link) 
	{
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


/**
 * Short one line description.                           (1)
 *
 * generates a set of site labels starting at n for gathering links
 * loop over even parity sites and gather with shift n from them 
 *
 * @param  n    -- 
 * @param  ind[] --
 */
    public static void makeindex(int n,int[] ind)
    {
	int[] x = new int[Global.DIM] ;

	split(x,n);

	int site=0;
	for(int iv=0;iv<vectorlength;iv++)
	    {
		while (parity[site] != 0 ) site++;

		ind[iv] = vshift(site,x);
		site++;
	    }
	return;
    }



    /**
     * Shift a site site_start by vector x[]
     *
     *
     * @param  site_start --  starting site
     * @param  x[] shift direction
     * @return site (integer)
     */
    public static int vshift(int site_start, int[] x)
    {
	int[] y = new int[Global.DIM];
	
	split(y,site_start);

	for(int i=0 ; i<Global.DIM ; i++)
	    {
		if (x[i] != 0)
		    {
			y[i] += x[i];

			while (y[i]>=Global.shape[i])
			    y[i]-=Global.shape[i];
			while (y[i]<0)
			    y[i]+= Global.shape[i];
		    }
	    }

	return siteindex(y);

}


   

    /**
     * Gives a unique index to site located at x[DIM] 
     *
     *
     * @param  x[] -- lattice coordinates
     * @return int -- site index
     */
    public static int siteindex(int[] x)
    {
	int result=0;

	for(int i=0;i<Global.DIM;i++)
	    result += shift[i]*x[i];

	return result;
    }



    /**
     * returns index of a site shifted dist in direction dir
     *
     * Longer description. If there were any, it would be    [2]
     * here.
     *
     * @param  site_start
     * @param  direction   -- direction that the sift is in
     * @param  dist
     * @return site (integer)
     */

    public static int ishift(int site_start, int dir, int dist)
    {
	int[] x = new int[Global.DIM];
	for(int i=0;i<Global.DIM;i++)
	    x[i]=0;

	x[dir]=dist;
	return vshift(site_start,x);
}




/**
 * Gather same color links into vector g starting at site 
 *
 * Longer description. If there were any, it would be    [2]
 * here.
 *
 * @param  gaugefield[] lattic
 * @param  int site -- ???????
 * @param  int link -- direction of gauge configuration
 * @return gaugefield[vectorlength]
 */
    public static gaugefield[]  getlinks(gaugefield[] lattice,int site,int gdir)
    {
	gaugefield[] g = new gaugefield[vectorlength] ;
	for(int i=0 ; i < vectorlength ; ++i)
	    {
		g[i] = new gaugefield(N);
	    }

	//	System.out.println("getlinks::site = " + site + " gdir = " + gdir);
	makeindex(site,myindex);

	int shift=nsites*gdir;
	for(int iv=0;iv<vectorlength;iv++)
	    g[iv]=lattice[myindex[iv]+shift].copy() ;
	

	return g ;
    }






/**
 * Short one line description.                           (1)
 *
 * accept new for old using metropolis algorithm
 * return average exponential of action change, this should fluctuate
 * about unity when in equilibrium
 * bias multiplies actions in exponential (i.e. beta/N)  
 * accepted changes returned in accepted
 *     actions passed in global variables sold and snew 
 *
 * @param  double bias        --
 * @param  gaugefield[] trial --
 * @param  gaugefield[] old   --
 * @return double 
 */

    public static double metro(gaugefield[] old,gaugefield[] trial,double bias) 
    {
	double expdeltas = 0.0 ;

	for(int iv=0;iv<vectorlength;iv++)
	    {
		double temp =  Math.exp((bias*(snew[iv]-sold[iv])));
		expdeltas += temp;
		
		if(generator.nextDouble()  < temp)
		    accepted[iv] = 1 ;
		else
		    accepted[iv] = 0 ;
		
		//		System.out.printf("DEBUG Accepted[%d] = %d\n",iv, accepted[iv]);
		//		System.out.printf("DEBUG tmp = %g Accepted[%d] = %d\n",temp, iv, accepted[iv]);

	    }

	for(int iv=0 ; iv<vectorlength ; iv++)
	    if (accepted[iv] != 0 ) 
		{
		    sold[iv] = snew[iv] ;
		    old[iv]  = trial[iv].copy()  ;
		}


  return expdeltas/vectorlength;
}   




    /**
     * scatter alternate links from vector g 
     *
     * Copy g into gauge configuration.
     *
     * @param  gaugefield[nlinks] lattice -- gauge configuration
     * @param  gaugefield[vectorlength] g 
     * @param  int site -- 
     * @param  int link -- link direction 
     */

    public static void savelinks(gaugefield[] lattice, gaugefield[]  g, int site, int link)
    {

	makeindex(site,myindex);
	int shift=nsites*link;
	
	for(int iv=0;iv<vectorlength;iv++)
	    {
		lattice[myindex[iv]+shift] = g[iv].copy();
	    }
	
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

	result=result/(N*vectorlength*count);

	System.out.printf("W[%d,%d] = %g\n",x,y,result);
	return result;
    }

    //
    //  DEBUG routines
    //


    public static void debug_init(int [] shift)
    {

	System.out.println("Dump out some tables\n");
	for (int i=0;i<shift.length;i++)
	    System.out.printf("Shift[%d] = %d\n" ,i, shift[i]   );

	for (int iv=0;iv<nsites;iv++)
	    {
		System.out.printf("Parity[%d] = %d\n" ,iv, parity[iv]   );
	    }

    }


    //
    // variables in the class
    //

/**
 * Description of the variable here.
 */
    public static int N ;

/**
 * Description of the variable here.
 */
    public static int nsites ;

/**
 * Description of the variable here.
 */
    public static int nlinks ; 

/**
 * Description of the variable here.
 */
    public static int nplaquettes ;

/**
 * Description of the variable here.
 */
    public static int vectorlength; 


    public static int[] accepted ;
    public static int[] myindex ;
    public static int[] parity ;

    public static double[] sold ;
    public static double[] snew ;

    public static int[] shift ;

/**
 * The gauge configurations as a one dimensional array
 */
    public static gaugefield[] ulinks; 

/**
 * A vector of random matrices
 */
    public static gaugefield[] table1 ;

/**
 * A vector of random matrices
 */
    public static gaugefield[] table2 ;

    public static gaugefield[] mtemp0 ;
    public static gaugefield[] mtemp1 ;
    public static gaugefield[] mtemp2 ;
    public static gaugefield[] mtemp3 ;
    public static gaugefield[] mtemp4 ;

/**
 * Random number generator 
 */
    private static Random generator ;

}
