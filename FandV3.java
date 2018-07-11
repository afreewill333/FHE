import java.io.*;
import java.util.*;
import java.math.*;

import org.mathIT.algebra.PolynomialZ;


public class FandV3{
	private static Random random = new Random();
	private static BigInteger ZERO = BigInteger.ZERO;
	private static BigInteger ONE = BigInteger.ONE;
	private static BigInteger TWO = new BigInteger("2");
		
	// sample a integer from [low, high]
	private static int randomInt(int low, int high)
	{   
        return (random.nextInt(high-low+1) + low);
    }
	
	/**
	 * BigInteger(int numBits, Random rnd) Constructs a randomly generated BigInteger, uniformly distributed over the range [0,   2^numBits - 1].
	 * using this function to get a BigInteger from Rq = {n: -q/2<n<=q/2 } = (-q/2, q/2] = [-q/2+1, q/2].
	 * @param The BigInteger q
	 * @return a BigInteger
	 */
	private static BigInteger RqSampling(BigInteger q)
	{
		///get the qn of q = 2^qn from q.///
		int qn;
		if(q.compareTo(ZERO)<0)
			qn = q.bitLength();
		else
			qn = q.bitLength() - 1;
		
		///step 1. get a random BigInteger bi from [0, 2^qn).///		
		BigInteger bi = new BigInteger(qn,random);
		///step 2. if bi > q/2 : bi -= q . that make bi shift to [-q/2+1, q/2].///
		if( bi.compareTo(q.divide(TWO)) > 0) 
			bi = bi.subtract(q);
		
		return bi;
	}
	
	/**
	 * @param mu(μ)  : a mean μ(mu);
	 * @param tau(τ) : a tail-cut parameter τ(tau)
	 * @param sigma  : a specific standard deviation σ(sigma)
	 * @return a integer from a bounded discrete Gaussian distribution draw.
	 */
	private static int mu = 0;
	private static int tau = 10;
	private static int GaussianSampling(int sigma)
	{	
		int Xmax = mu + tau*sigma;
		int Xmin = mu - tau*sigma;

		while(true)
		{
			int x = randomInt(Xmin,Xmax);
			double p = 1.0/Math.sqrt(2*Math.PI)/sigma * Math.pow(Math.E, -0.5*(x-mu)/sigma*(x-mu)/sigma);
			double r = random.nextDouble();
			if(r < p) return x;
		}

	}
	
	/**
	 * @param res  : the polynomial
	 * @param q    : from Rq = Zq[x]/... . 
	 * @return the res with coefficients that have been Centered. 
	 */
	private static PolynomialZ CentredCoeff(PolynomialZ res, BigInteger q)
	{
		/// res is a map of < exponent, coefficient >. ///
		for(BigInteger idx : res.keySet())
		{
			BigInteger value = res.get(idx);
			value = value.mod(q);
			
			/// get i = q/2. ///
			BigInteger i;// i=n for q=2n+1or2n   BigInteger.getLowestSetBit() 
			if(q.mod(BigInteger.TWO).equals(BigInteger.ZERO))
				i = q.divide(BigInteger.TWO); // or just BigInteger i = q.divide(BigInteger.TWO);
			else
				i = q.subtract(BigInteger.ONE).divide(BigInteger.TWO);
			
			/// if value>q/2 : value -= q. so that the coefficient shifts. ///
			if(value.compareTo(i) > 0)
				value = value.subtract(q);
			
			res.put( idx, value );
		}
		
		return res;
	}
	
	private static PolynomialZ MOD(PolynomialZ x, PolynomialZ y)
	{
		BigInteger degX = x.getDegree();
		BigInteger degY = y.getDegree();
		
		if(degX.compareTo(degY)>0)
			return x.mod(y);
		else
		if(degX.compareTo(degY)==0)
		{
			if(x.get(degX).compareTo(y.get(degY))>=0) // also need: [x.get(degX)] %  [y.get(degY)] == 0
				return x.mod(y);
			else
				return x; /// will not happen in the case of g(x)%(x^n + 1)!
			}
		else
		if(degX.compareTo(degY)<0)
			return x;
		
		return new PolynomialZ();
	}
	private static void print(PolynomialZ r) {
		
		//System.out.println("- - - - - - - - - - - - - - - - - - - - - -");
		int c=0;
		System.out.print("( ");
		String gap = "    ";
		List<BigInteger> idxs = new ArrayList<BigInteger>(r.keySet());
		for(BigInteger i : idxs )
		{
			if(r.get(i)!=null && r.get(i).compareTo(ZERO)!=0)
			{
				if(i.equals(ZERO)) {
					System.out.print(r.get(i));
					continue;
				}
				if(i.equals(ONE)) {
					if(r.get(i).compareTo(ONE)!=0)
					{
						System.out.print(r.get(i)+ "x" + gap + "+" + gap );
					    continue;
					}else {
						System.out.print("x" + gap + "+" + gap );
					    continue;
					}
				}
				
				if(r.get(i).equals(ONE)) 
					System.out.print("x^" + i + gap + "+" + gap);
				else   
					System.out.print(r.get(i)+ "x^" + i + gap + "+" + gap);
				
				c++;
				if(c>12) break;
			}
		}
		System.out.print(" )");
		//System.out.println("- - - - - - - - - - - - - - - - - - - - - -");
	}

	public static void main(String[] args){
		// 1s = 10^3ms = 10^6us = 10^9ns
		long beg = System.nanoTime();

		/**
		 * Notation:
		 *   q     : modulus in the cipher text space(coefficient modulus) of the form q1*q2*...*qk,where qi are prime
		 *   t     : modulus in the plain text space(plain text modulus)
		 *   x^n+1 : The polynomial modulus which specifies the ring R
		 *   R     : The ring Z[x]/f(x), where f(x) =(x^n+1)=x^(2^(d-1))+1
		 *   Ra    : The ring Za[x]/(x^n+1),i.e. same as the ring R but with coefficients reduced modulo a
		 *   det   : Quotient on division of q by t, or [q/t]
		 *   rt(q) : Remainder on division of q by t, i.e. q = det*t + rt(q), where 0<=rt(q)<t
		 *   X     : Error distribution(a truncated discrete Gaussian distribution)
		 *   sigma : Standard deviation of X
		 *   B     : Bound on the distribution X
		 *   
		 *   The scheme FV : SecretKeyGen,PublicKeyGen,EvaluationKeyGen,Encrypt,Decrypt,Add and Multiply.
		 */
		int d = 13;//4095 degree polynomials)
		BigInteger q = TWO.pow(128);
		BigInteger t = TWO.pow(15);
		int sigma = 16;
		PolynomialZ Phi = new PolynomialZ();
		Phi.put(TWO.pow(d-1), ONE);
		Phi.put(ZERO, ONE);
		System.out.println("Phi is done.");
		///SecretKeyGen(lambda): Sample s <- R2 and output sk = s.
		PolynomialZ sk = new PolynomialZ();
		long length = (long) Math.pow(2, d-1);
		for(long i=0;i<length;i++) {
			int r = randomInt(0,1);
			sk.put(BigInteger.valueOf(i), BigInteger.valueOf(r));
		}
		PolynomialZ s = sk;
		System.out.println("sk is done.");
		///PublicKeyGen(sk): Set s = sk, sample a <- Rq, and e <- X. Output pk = ([-(a*s+e)]q, a).
		PolynomialZ a = new PolynomialZ();
		for(long i=0;i<length;i++) {
			BigInteger coefficient = RqSampling(q);
			a.put(BigInteger.valueOf(i), coefficient);
		}
		PolynomialZ e = new PolynomialZ();
		for(long i=0;i<length;i++) {
			BigInteger coefficient = BigInteger.valueOf(GaussianSampling(sigma));
			e.put(BigInteger.valueOf(i), coefficient);
		}
		
		PolynomialZ p0 = a.multiply(s); MOD(p0, Phi);
		p0 = p0.plus(e);PolynomialZ zero = new PolynomialZ();
		p0 = zero.minus(p0);
		CentredCoeff(p0, q);
		PolynomialZ p1 = a;
		System.out.println("pk is done.");
		///EvaluationKeyGen(sk,w): For i{0,...,l},sample ai <- Rq, ei <- X. Output
		///                        evk = ([-(ai*s + ei) + w^i*s^2]q, ai).
		
		///Encrypt(pk,m): For m in Rt, let pk=(p0,p1).Sample u <- R2, and e1,e2 <- X. Compute
		///                        ct = ([det*m + p0*u +e1]q, [p1*u + e2]q)
		BigInteger message =  BigInteger.valueOf(3636);
		PolynomialZ mx = new PolynomialZ();
		long deg = 0;
		for(char ch : message.toString(2).toCharArray())
		{////////////////should rewrite!///////////////////
			mx.put(BigInteger.valueOf(message.toString(2).length()-deg-1), BigInteger.valueOf(ch-'0'));// work when base<10
			deg++;
		}
		PolynomialZ det = new PolynomialZ(ZERO, q.divide(t));
		PolynomialZ u  = new PolynomialZ();
		PolynomialZ e1 = new PolynomialZ();
		PolynomialZ e2 = new PolynomialZ();
		for(long i=0;i<length;i++) {
			u.put(BigInteger.valueOf(i), BigInteger.valueOf(randomInt(0,1)));
			e1.put(BigInteger.valueOf(i), BigInteger.valueOf(GaussianSampling(sigma)));
			e2.put(BigInteger.valueOf(i), BigInteger.valueOf(GaussianSampling(sigma)));
		}
		
		PolynomialZ c0 = det.multiply(mx);
		PolynomialZ p0u = p0.multiply(u);
		p0u = MOD(p0u, Phi);
		c0 = c0.plus(p0u).plus(e1);
		CentredCoeff(c0, q);
		System.out.println("c0 is done.");
		PolynomialZ c1 = p1.multiply(u);
		c1 = MOD(c1, Phi);
		c1 =c1.plus(e2);
		CentredCoeff(c1, q);
		System.out.println("c1 is done.");
		
		///Decrypt(sk,ct): Set s = sk, c0 = ct[0], and c1 = ct[1]. Output
		PolynomialZ result = c1.multiply(s);
		result = MOD(result, Phi);
		result = c0.plus(result);
		CentredCoeff(result, q);
		
		BigDecimal tlq = new BigDecimal(t).divide(new BigDecimal(q));
		for(BigInteger idx : result.keySet()) {
			BigDecimal coefficient = tlq.multiply( new BigDecimal(result.get(idx)) );
			coefficient = coefficient.setScale(0, RoundingMode.HALF_DOWN);
			result.put(idx, coefficient.toBigInteger());			
		}
		CentredCoeff(result, t);
		
		System.out.println("m(2) = " + result.evaluate(TWO));
		

		
		
		
		
		
		
		long end = System.nanoTime();
		System.out.println();
		System.out.println("end-beg = " + (end-beg)/1e6 + "ms = " + (end-beg)/1e9 +"s." );


	}

}