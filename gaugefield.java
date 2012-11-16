/**
 * @author      Craig McNeile
 * @version     0.5      
 * @since       2012-11-16  
 */

public class gaugefield 
{

    public gaugefield (int d)
    {
	this.N = d ;
	real = new double[d][d] ;
	imag = new double[d][d] ;
    }

    public void printmatrix()
    {
	int i,j;
	for (i=0;i<N;i++){
	    System.out.printf("\n");
	    for (j=0;j<N;j++){
		System.out.printf(" (%g, %g)   ",real[i][j],imag[i][j]);
	    }
	}
	System.out.printf("\n");
    }


    /* projects a matrix onto the group SU(N) */
    public void project()
    {
	int i,j,k,nmax;
	//	nmax=N-(N<4); /* 2 and 3 are treated specially */

	if( N<4 )
	    nmax = N - 1  ;
	else
	    nmax = N ;

	/* loop over rows */
	for (i=0;i<nmax;i++) {
	    /* normalize i'th row */
	    double temp=real[i][0] * real[i][0]
		+ imag[i][0]*imag[i][0];
	    for (j=1;j<N;j++)
		temp+= real[i][j]*real[i][j]
		    +imag[i][j]*imag[i][j];
	    temp=1/Math.sqrt(temp);   
	    for (j=0;j<N;j++) {
		real[i][j]*=temp;
		imag[i][j]*=temp;
	    }
	    /* orthogonalize remaining rows */
	    double adotbr,adotbi;
	    for (k=i+1;k<nmax;k++) {
      adotbr=real[i][0]*real[k][0]
	+imag[i][0]*imag[k][0];
      adotbi= real[i][0]*imag[k][0]
	-imag[i][0]*real[k][0];
      for (j=1;j<N;j++) {
	adotbr+=real[i][j]*real[k][j]
	  + imag[i][j]* imag[k][j];
	adotbi+= real[i][j]*imag[k][j]
	  -imag[i][j]*real[k][j];
      }
      for (j=0;j<N;j++) {
	real[k][j]-=adotbr* real[i][j]
	  -adotbi*imag[i][j];
	imag[k][j]-= adotbr*imag[i][j]
	  +adotbi* real[i][j];
      } 
    } /* end of k loop */
  } /* end of i loop */
    /* remove determinant, treating group=2 or 3 as special cases */

  switch (N) {
  case 3:
    thirdrow();
    break; 
  case 2: /* for su(2) */
    real[1][0]= -real[0][1];
    real[1][1]=  real[0][0];
    imag[1][0]=  imag[0][1];
    imag[1][1]= -imag[0][0];
    break;  
  default: /* remove the determinant from the first row */
    double x=0.0 , y=0.0 ,w;
    this.determinant(x,y);
    for (i=0;i<N;i++) {
      w=real[0][i]*x
	+ imag[0][i]*y;
      imag[0][i]= imag[0][i]*x
	- real[0][i]*y;
      real[0][i]=w;
    }
  } /* end switch */

    }



    public gaugefield  conjugate() 
    {
	gaugefield result = new gaugefield(N) ;

	for(int i=0;i<N;i++)
	    for(int j=0;j<N;j++) 
		{
		    result.real[i][j]= real[j][i];
		    result.imag[i][j]= -imag[j][i];
		}


	return result;
    }


    public void set_constant(double re, double im)
    {
	int i,j;
	for (i=0;i<N;i++)
	    {
		for (j=0;j<N;j++)
		    {
			real[i][j] = re ; 
			imag[i][j] = im ;
		    }
	    }


    }



    public void set_unit()
    {
	int i,j;
	for (i=0;i<N;i++){
	    for (j=0;j<N;j++){
		real[i][j] = 0.0 ; 
		imag[i][j] = 0.0 ;

	    }
	}

	for (i=0;i<N;i++)
	    real[i][i] = 1.0 ; 

    }



    public double trace_re()
    {
	double ans = 0.0  ;

	int i ;
	for (i=0;i<N;i++)
	    ans += imag[i][i]  ;

	return ans ;

    }


    public gaugefield copy () {
	gaugefield X = new gaugefield(N);

	for (int i = 0; i < N ; i++) {
	    for (int j = 0; j < N  ; j++) {
		X.real[i][j] = real[i][j] ;
		X.imag[i][j] = imag[i][j] ;
	    }
	}
	return X;
    }




    public gaugefield prod (gaugefield Y) 
    {
	gaugefield X = new gaugefield(N);

	for (int i = 0; i < N ; i++) {
	    for (int j = 0; j < N  ; j++) {
		X.real[i][j] = 0.0 ;
		X.imag[i][j] = 0.0 ;

	    for (int k = 0; k < N  ; k++) {
		X.real[i][j] += real[i][k] * Y.real[k][j] ;
		X.real[i][j] -= imag[i][k] * Y.imag[k][j] ;

		X.imag[i][j] += real[i][k] * Y.imag[k][j]  ;
		X.imag[i][j] += imag[i][k] * Y.real[k][j]  ;
	    }


	    }
	}
	return X;
    }



/**
 * Check the Unitarity of the gauge matrix
 *
 * Check that U^dagger U = 1 
 *
 *
 */
    public void check_unitarity () 
    {
	for (int i = 0; i < N ; i++) 
	    {
		for (int j = 0; j < N  ; j++) 
		    {
			double Tre = 0.0 ;
			double Tim = 0.0 ;
			for (int k = 0; k < N  ; k++) 
			    {
				Tre += real[i][k] * real[j][k] ;
				Tre += imag[i][k] * imag[j][k] ;

				Tim -= real[i][k] * imag[k][j]  ;
				Tim += imag[i][k] * real[k][j]  ;
			    }
			System.out.printf("UU^dagger[%d,%d] = (%g, %g)   ",i,j, Tim, Tim);
	    }
	}
    }





    public gaugefield sub (gaugefield Y) 
    {
	gaugefield X = new gaugefield(N);

	for (int i = 0; i < N ; i++) {
	    for (int j = 0; j < N  ; j++) {
		X.real[i][j] = real[i][j] - Y.real[i][j] ;
		X.imag[i][j] = imag[i][j] - Y.imag[i][j]  ;
	    }
	}
	return X;
    }




    //    need to return 
    //    public void determinant(double& detr, double& deti) {

    public void determinant(double detr, double deti) {
	/* subroutine to calculate matrix determinant */
	/* could perhaps improve later by doing group=2 and 3 by hand,
	   but as of now these cases don't call this routine anyway */
	int i,j,k,im,km;
	double temp;
	
	im=N-1;
	for (i=0;i<im;i++) {
	    km=i+1;
	    /* find magnitude of i'th diagonal element */
	    temp=1./(this.real[i][i]*this.real[i][i]
	     +this.imag[i][i]*this.imag[i][i]);
	    for (k=km;k<N;k++) {
		/* subtract part of row i from row k to make [k][i] element vanish */
		/* inner product of the rows */
		detr=(this.real[k][i]*this.real[i][i]
		      +this.imag[k][i]*this.imag[i][i])*temp;
		deti=(this.imag[k][i]*this.real[i][i]
		      -this.real[k][i]*this.imag[i][i])*temp;
		for (j=km;j<N;j++) {
		    this.real[k][j]+= -detr*this.real[i][j]
			+deti*this.imag[i][j];
		    this.imag[k][j]+=  -detr*this.imag[i][j]
			-deti*this.real[i][j];
		}
	    }  
	}  
	/* multiply diagonal elements */
	detr=this.real[0][0];
	deti=this.imag[0][0];
	for (i=1;i<N;i++) {
	    temp=detr*this.real[i][i]-deti*this.imag[i][i];
	    deti=deti*this.real[i][i]+detr*this.imag[i][i];
	    detr=temp;
	}
	return;
    }

    private void thirdrow() 
    {
	/* for su(3) construct third row from first two */
	int i,j,k;
	for (i=0;i<3;i++) {
	    j=(i+1)%3;       
	    k=(i+2)%3;
	    real[2][i]= real[0][j]*real[1][k]
		-imag[0][j]*imag[1][k]
		-real[1][j]*real[0][k]
		+imag[1][j]*imag[0][k];
	    imag[2][i]=-real[0][j]*imag[1][k]
		-imag[0][j]*real[1][k]
		+real[1][j]*imag[0][k]
		+imag[1][j]*real[0][k];
	}

}


    //    gaugefield& project(); /* projects onto the gauge group */
    // void determinant(double& rez, double& imz);
    // void printgaugefield();
    // gaugefield& operator= (double x);
    // gaugefield& operator*= (double x);


/**
 * Real part of SU(N) matrix
 */
    public double[][] real ;

/**
 * Imaginary part of SU(N) matrix
 */
    public double[][] imag ;

/**
 * Group SU(N)
 */
    private int N ;

};

