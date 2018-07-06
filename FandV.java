import java.io.*;
import java.util.*;
import java.math.*;

import org.mathIT.algebra.PolynomialZ;

public class FandV{
	private static Random random = new Random();
	private static BigInteger TWO = new BigInteger("2");
	private static BigInteger ONE = BigInteger.ONE;
	private static BigInteger ZERO = BigInteger.ZERO;

	// sample a integer from [low, high]
	private static int randomInt(int low, int high)
	{   
        return (random.nextInt(high-low+1) + low);
    }
 
	// sample a BigInteger from [ -(2^numBits-1), (2^numBits-1)] ; 
	// return a BigInteger; 
	private static BigInteger RqSampling(int numBits)
	{
		BigInteger bi = new BigInteger(numBits,random);
		if(random.nextBoolean()) 
			return bi;
		else
			return bi.negate();
	}
	// a mean μ(mu); 
	// a specific standard deviation σ(sigma); 
	// a tailcut parameter τ(tau); 
	private static int mu = 0;
	private static int tau = 10;
	private static int GaussianSampling(int sigma)
	{	
		double h = -Math.PI/sigma/sigma;
		int Xmax = mu + tau*sigma;
		int Xmin = mu - tau*sigma;

		while(true)
		{
			int x = randomInt(Xmin,Xmax);
			double p = Math.pow(Math.E, h*(x-mu)*(x-mu));
			double r = random.nextDouble();
			if(r < p) return x;
		}

	}
	//*******************************************
	//*  module the coeff of Zq[x] where q=17.  *
	//*******************************************
	private static PolynomialZ CentredCoeff(PolynomialZ res, BigInteger bp)
	{
		//BigInteger bp = new BigInteger("17");
		//System.out.println("~ ~ ~ ~ ~ ~ ~ ~");
		for(BigInteger bi:res.keySet())
		{
			BigInteger v = res.get(bi).mod(bp);
			BigInteger i;// i=n for bp=2n+1or2n
			if(bp.mod(BigInteger.TWO).equals(BigInteger.ZERO))
				i = bp.divide(BigInteger.TWO);
			else
				i = bp.subtract(BigInteger.ONE).divide(BigInteger.TWO);
			if(v.compareTo(i)>0)
				v = v.subtract(bp);
			res.put( bi, v );
		}
		//!!!!!!!!!!!!!!!!!!!!
		return res;
	}


	public static void main(String[] args){
		// 1s = 10^3ms = 10^6us = 10^9ns
		long beg = System.nanoTime();

		// Notation:
		//   Zq[x] denote polynomial whose coefficients belong to Zq.
		//   Zq is the set of integers {n: -q/2<n<=q/2 }.
		//   Rq[x] = Zq[x]/(x^(2^(d-1)) + 1)
		//   ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
		//   d = 13 (4095 degree polynomials)
		//   q = 2^128
		//   t = 2^15
		//   σ = 16 (sigma)
		int d = 13;
		BigInteger q = TWO.pow(128); int qn = 128;
		BigInteger t = TWO.pow(15); int tn = 15;
		int sigma = 16;
		//Stage.1 Key Generation
		//  Step.1 The secret key, Ks, is simply a uniform random draw from R2(-1,1].
		//         sample a 2^(d-1) binary vector for the polynomial coefficients.
		PolynomialZ Ks = new PolynomialZ();
		//Ks.put( n, ai );	fx = ai*x^n	 
		for(BigInteger i=ZERO;i.compareTo(TWO.pow(d-1))<0;i=i.add(ONE)){
			Ks.put(i, BigInteger.valueOf(random.nextInt(2)));
			//System.out.println(Ks.get(i) + " x^ " + i);
		}
		//  Step.2 The public key, Kp, is a vector containing two polynomials:
		//         Kp = (Kp1, Kp2) = ([-(a*Ks + e)]%q, a)   a ~ Rq, e ~ X(0,σ)
		PolynomialZ a = new PolynomialZ();	 
		for(BigInteger i=ZERO;i.compareTo(TWO.pow(d-1))<0;i=i.add(ONE)){
			a.put(i, RqSampling(qn-1));
		}
		PolynomialZ e = new PolynomialZ();	 
		for(BigInteger i=ZERO;i.compareTo(TWO.pow(d-1))<0;i=i.add(ONE)){
			e.put(i, BigInteger.valueOf(GaussianSampling(sigma)));
		}
		PolynomialZ fx = new PolynomialZ();	 
		fx.put(TWO.pow(d-1), BigInteger.ONE);
		fx.put(BigInteger.ZERO, BigInteger.ONE);
	
		PolynomialZ Kp1;
		PolynomialZ neg = new PolynomialZ(ZERO, ONE.negate());
		Kp1 = a.multiply(Ks).plus(e).multiply(neg);		
		Kp1 = Kp1.mod(fx);
		CentredCoeff(Kp1,BigInteger.valueOf(qn-1));
		PolynomialZ Kp2 = a;//(PolynomialZ) a.clone();	 

		//Stage.2 Encryption
		//  Step.1 An integer message m is first represented as m(x)∈ Rt.
		BigInteger m = new BigInteger("210281198811199375");
		String mxs = m.toString(2);
		PolynomialZ mx = new PolynomialZ();	 
		for(int i=0;i<mxs.length();i++){
			mx.put(BigInteger.valueOf(i), new BigInteger(""+(mxs.charAt(mxs.length()-i-1))) );
		}
		//System.out.println("m = " + mx.evaluate(TWO) );
		
		//  Step.2 Encryption then renders a cipher text which is a vector containing two polynomials:
		//         c = (c1,c2) = ( [Kp1*u + e1 + det*m(x)]q, [Kp2*u + e2]q )
		//         u,e1,e2 ~ X(0,σ)   det = [q/t] 
		PolynomialZ det = new PolynomialZ(BigInteger.ZERO, q.divide(t));
		PolynomialZ u = new PolynomialZ();
		PolynomialZ e1= new PolynomialZ();
		PolynomialZ e2= new PolynomialZ();
		for(BigInteger i=ZERO;i.compareTo(TWO.pow(d-1))<0; i=i.add(ONE)){
			u.put(i, BigInteger.valueOf(GaussianSampling(sigma)));
			e1.put(i, BigInteger.valueOf(GaussianSampling(sigma)));
			e2.put(i, BigInteger.valueOf(GaussianSampling(sigma)));
		}
		PolynomialZ c1 = Kp1.multiply(u).plus(e1).plus(det.multiply(mx));
		c1 = c1.mod(fx); CentredCoeff(c1, q);
		PolynomialZ c2 = Kp2.multiply(u).plus(e2);
		c2 = c2.mod(fx); CentredCoeff(c2, q);

		//Stage.3 Decryption
		//  Step.1 decryption of a cipher text c is by evaluating:
		//              m0(x) = [ [t*[c1 + c2*Ks]q / q ] ]t
		//         so that m = m0(x).				
		PolynomialZ m0 = c1.plus(c2.multiply(Ks)).mod(fx);
		for(BigInteger i=ZERO;i.compareTo(m0.getDegree())==0; i=i.add(ONE))
		{
			BigInteger coef = m0.get(i).multiply(t).divide(q);
			m0.put(i, coef);
		}
		CentredCoeff(m0, t);
		
		System.out.println("m0(2) = " + m0.evaluate(TWO));
		









		
		
		
		
		







		long end = System.nanoTime();
		System.out.println();
		System.out.println("end-beg = " + (end-beg)/1e6 + "ms = " + (end-beg)/1e9 +"s." );

		
		/*
		 * Test The Gaussian Sampling 
		Map<Integer,Integer> map = new TreeMap<Integer,Integer>();
		for(int i=0;i<1000;i++)
		{
			int kx = GaussianSampling(6);
			if(map.containsKey(kx))
				map.put(kx,map.get(kx)+1);
			else
				map.put(kx,1);
		}
		for(int k : map.keySet())
		{
			String bar = new String();
			for(int i=0;i<map.get(k)/10;i++) bar += "*";
			System.out.println((0<=k&&k<10?" ":"")+k + " : " + bar);
		}		
		*/

		/*
		 * 
		BigDecimal bd = new BigDecimal(Math.E);
		System.out.println("e^(-50) = " + BigDecimal.ONE.divide(bd.pow(50),128, RoundingMode.DOWN));
		System.out.println();
		System.out.println("Hello, World! "+e);
		for(int i=0;i<100;i++)
		{
			System.out.println(RqSampling(127));
		}
		*/
	}
}