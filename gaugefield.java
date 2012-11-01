public class gaugefield {


    public gaugefield (int d)
    {
	this.GROUP = d ;
	real = new double[d][d] ;
	imag = new double[d][d] ;
    }



    public void printmatrix(){
	int i,j;
	for (i=0;i<GROUP;i++){
	    System.out.printf("\n");
	    for (j=0;j<GROUP;j++){
		System.out.printf(" (%g, %g)   ",real[i][j],imag[i][j]);
	    }
	}
	System.out.printf("\n");
    }


    public void project()
    {




    }



    public gaugefield  conjugate() 
    {
	gaugefield result = new gaugefield(GROUP) ;

	for(int i=0;i<GROUP;i++)
	    for(int j=0;j<GROUP;j++) 
		{
		    result.real[i][j]= real[j][i];
		    result.imag[i][j]= -imag[j][i];
		}


	return result;
    }


    public void set_constant(double re, double im)
    {
	int i,j;
	for (i=0;i<GROUP;i++){
	    for (j=0;j<GROUP;j++){
		real[i][j] = re ; 
		imag[i][j] = im ;

	    }
	}


    }



    public double trace_re()
    {
	double ans = 0.0  ;

	int i ;
	for (i=0;i<GROUP;i++)
	    ans += imag[i][i]  ;

	return ans ;

    }


    public gaugefield copy () {
	gaugefield X = new gaugefield(GROUP);

	for (int i = 0; i < GROUP ; i++) {
	    for (int j = 0; j < GROUP  ; j++) {
		X.real[i][j] = real[i][j] ;
		X.imag[i][j] = imag[i][j] ;
	    }
	}
	return X;
    }




    public gaugefield prod (gaugefield Y) {
	gaugefield X = new gaugefield(GROUP);

	for (int i = 0; i < GROUP ; i++) {
	    for (int j = 0; j < GROUP  ; j++) {
		X.real[i][j] = 0.0 ;
		X.imag[i][j] = 0.0 ;

	    for (int k = 0; k < GROUP  ; k++) {
		X.real[i][j] += real[i][k] * Y.real[k][j] ;
		X.real[i][j] -= imag[i][k] * Y.imag[k][j] ;

		X.imag[i][j] += real[i][k] * Y.imag[k][j]  ;
		X.imag[i][j] += imag[i][k] * Y.real[k][j]  ;
	    }


	    }
	}
	return X;
    }





    public gaugefield sub (gaugefield Y) {
	gaugefield X = new gaugefield(GROUP);

	for (int i = 0; i < GROUP ; i++) {
	    for (int j = 0; j < GROUP  ; j++) {
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
	
	im=GROUP-1;
	for (i=0;i<im;i++) {
	    km=i+1;
	    /* find magnitude of i'th diagonal element */
	    temp=1./(this.real[i][i]*this.real[i][i]
	     +this.imag[i][i]*this.imag[i][i]);
	    for (k=km;k<GROUP;k++) {
		/* subtract part of row i from row k to make [k][i] element vanish */
		/* inner product of the rows */
		detr=(this.real[k][i]*this.real[i][i]
		      +this.imag[k][i]*this.imag[i][i])*temp;
		deti=(this.imag[k][i]*this.real[i][i]
		      -this.real[k][i]*this.imag[i][i])*temp;
		for (j=km;j<GROUP;j++) {
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
	for (i=1;i<GROUP;i++) {
	    temp=detr*this.real[i][i]-deti*this.imag[i][i];
	    deti=deti*this.real[i][i]+detr*this.imag[i][i];
	    detr=temp;
	}
	return;
    }



    //    gaugefield& project(); /* projects onto the gauge group */
    // void determinant(double& rez, double& imz);
    // void printgaugefield();
    // gaugefield& operator= (double x);
    // gaugefield& operator*= (double x);



    public double[][] real ;
    public double[][] imag ;

    private int GROUP ;

};

