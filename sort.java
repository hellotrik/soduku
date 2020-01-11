package ran.good;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
/*
 * @Sort
 */
public final class Sort {
	public interface Mable{
		void move(int n,int a,int b);
	}
	/**
	 * n!
	 */
	public static BigInteger npow(int n){
		BigInteger x=new BigInteger("1");
		if(n<2)return x;
		
		for(int i=2;i<=n;i++){
			x=x.multiply(new BigInteger(i+""));
		}
		return x;
	}
	/**
	 * steling 根号(2πn)*e^(theta/12n)*(n/e)^n
	 * @param n
	 * @return
	 */	
	public static BigInteger steling(int n){
		BigInteger x=new BigInteger("1");
		if(n<2)return x;
		double ret=0;
		ret=Math.sqrt(Math.PI*2*n);
		ret*=Math.exp(Math.random()/(12*n));
		BigDecimal y=new BigDecimal(n/Math.E+"");
		y=y.pow(n);
		y=y.multiply(new BigDecimal(ret+""));
		return y.toBigInteger();
	}
	
	//---------Einstan
	private static String[] st={"挪威","英国","瑞典","丹麦","德国","蓝","红","绿","白","黄","狗","鸟","猫","马","鱼",
			"Pr","BM","DH","PM","BL","牛奶","茶","啤酒","矿泉水","咖啡"};
	private static int[] t={0,3,4,5,6,1,3,9,10,11,4,12,13,14,15,6,7,11,12,16,2,5,7,8,9};
	private static void d(int[] x){
		for(int i=1;i<6;i++){
			System.out.print(i+"号房子: ");for(int j=0;j<25;j++)if(x[t[j]]==i) System.out.print(st[j]+",");System.out.println();
		}
	}
	private static boolean c(int[] x){
		if(x[9]>0&&x[10]>0&&x[9]+1!=x[10]) return false;
		if(x[11]>0&&x[14]>0&&Math.abs(x[11]-x[14])!=1) return false;
		if(x[16]>0&&x[13]>0&&Math.abs(x[16]-x[13])!=1) return false;
		if(x[16]>0&&x[8]>0&&Math.abs(x[16]-x[8])!=1) return false;
		for(int i=0;i<5;i++)for(int j=0;j<5;j++)if(x[t[i*5+j]]>0)for(int k=j+1;k<5;k++)if(x[t[i*5+j]]==x[t[i*5+k]])
				return false;
		return true;
	}
	public static void Einstan(boolean b){
		int[] x=new int[17]; for(int i=0;i<3;i++)x[i]=i+1;
		if(b) Einstan(x,3);	else Einstan(x);
	}
	private static boolean Einstan(int[] x,int dep){
		if(dep>16){d(x);return true;}
		for(int i=1;i<6;i++){
			x[dep]=i;
			if(c(x)) {if(Einstan(x,dep+1))return true;}
			x[dep]=0;
		}
		return false;
	}
	private static void Einstan(int[] x){
		int dep=3;	int[] f=new int[14];	for(int i=0;i<14;i++)f[i]=1;
		while(dep>2&&dep<17){
			if(f[dep-3]<6){
				for(int i=f[dep-3];i<6;i++){
					x[dep]=i;f[dep-3]=i+1;
					if(c(x)&&dep<17){
						if(dep==16){d(x);return;}
						else {dep++;break;}
					}
				}
			}else{
				f[dep-3]=1;x[dep]=0;dep--;
			}
		}
	}
	
	//----bsearch
	public static <T extends Comparable<? super T>> int bsearch(T[] a,T x){
		int lf=0,r=a.length-1;
		while(lf<=r){
			int mid=(lf+r)>>>1;
			if(x==a[mid]) return mid;
			if(x.compareTo(a[mid])==1)lf=mid+1;
			else r=mid-1;
		}
		return -1;
	}
	//---------HanoiTower
	public static void Hanoi(int n){
		Mable m=new Mable(){
			String[] s={"A","B","C"};
			@Override
			public void move(int n, int a, int b) {
				System.out.println("Tower"+s[a]+"-"+n+"--->Tower"+s[b]);
				
			}
		};
		Hanoi(n, 0,1,2,m);
	}
	public static void Hanoi(int n,int a,int b,int c,Mable moveable){
		if(n>0){
			Hanoi(n-1,a,c,b,moveable);
			moveable.move(n, a, b);
			Hanoi(n-1,c,b,a,moveable);
		}
	}
	/*Z+划分
	 * 整数n的所有划分中,最大加数不大于m的划分记为q(n,m)
	 */
	public static int q(int n){return q(n,n);}
	private static int q(int n,int m){
		if(n<1||m<1) return 0;
		if(n==1||m==1)return 1;
		if(n<m) return q(n,n);
		if(n==m)return q(n,m-1)+1;
		return q(n,m-1)+q(n-m,m);
	} 
	//--------排列
	public static <T extends Comparable<? super T>>
	void perm(ArrayList<T[]> res,T[] a,int k,int m){
		if(k==m)res.add(a.clone());
		else for(int i=k;i<=m;i++){
			swap(a,k,i);
			perm(res,a,k+1,m);
			swap(a,k,i);
		}
	}
	//------读取字符串
	public static String gets(){
		BufferedReader br=new BufferedReader(new InputStreamReader(System.in));
		String s="";
		try {
			s = br.readLine();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return s;
	}
	//-------------fb数
	/**
	 * <h>不建议使用此线性复杂度函数
	 * 建议使用{@link #Fib(long)}对数复杂度函数</h>
	 * 
	 */
	public static BigInteger Fibold(long n){
		if(n<=1) return new BigInteger("1");
		BigInteger sn=new BigInteger("1"),sn_1=new BigInteger("0");
		for(long i=0;i<n;i++){
			BigInteger t=sn;
			sn=sn.add(sn_1);
			sn_1=t;
		}
		return sn;
	}
	
	public static BigInteger Fib(long n){
		if(n<2)return new BigInteger("1");
		return Mpow(n)[0];
	}
	
	public static BigInteger[] Mpow(long n){
		BigInteger[] m={new BigInteger("1"),new BigInteger("1"),new BigInteger("1"),new BigInteger("0")};
		BigInteger[] m1={new BigInteger("1"),new BigInteger("1"),new BigInteger("1"),new BigInteger("0")};
		if(n==1) return m;
		if((n&1)==0){m=Mpow(n>>>1);m=M2x2(m,m);}
		else {
			m=Mpow((n-1)>>>1);
			m=M2x2(m,m);
			m=M2x2(m,m1);
		}
		return m;
	}
	
	public static BigInteger[] M2x2(BigInteger[] m,BigInteger[] m2){
		BigInteger[] c=new BigInteger[4];
		c[0]=m[0].multiply(m2[0]).add(m[1].multiply(m2[2]));
		c[1]=m[0].multiply(m2[1]).add(m[1].multiply(m2[3]));
		c[2]=m[2].multiply(m2[0]).add(m[3].multiply(m2[2]));
		c[3]=m[2].multiply(m2[1]).add(m[3].multiply(m2[3]));
		return c;
	}
	
	
	//------radixsort
	public static void radixSortA(String[] a,int len){
		final int BKTS=256;
		ArrayList<ArrayList<String>> bs=new ArrayList<>();
		for(int i=0;i<BKTS;i++)bs.add(new ArrayList<String>());
		for(int p=len-1;p>=0;p--){
			for(String s:a)bs.get(s.charAt(p)).add(s);
			int idx=0;
			for(ArrayList<String> tb:bs){
				for(String s:tb)
					a[idx++]=s;
				tb.clear();
			}
		}
	}
	
	//--------heapsort
	private static int lc(int i){return i<<1;}
	private static <T extends Comparable<? super T>>
	void downheap(T[] a,int i,int n){
		int c;T t;
		for(t=a[i];lc(i)<n;i=c){
			c=lc(i);
			if(c!=n&&a[c].compareTo(a[c+1])<0)c++;
			if(t.compareTo(a[c])<0)a[i]=a[c];
			else break;
		}
		a[i]=t;
	}
	public static <T extends Comparable<? super T>>
	void heapSort(T[] a){
		 for (int i = (a.length>>>1)-1; i >=0; i--)
	            downheap(a,i,a.length-1);
		 for(int i=a.length-1;i>0;i--){swap(a,0,i);	downheap(a,0,i-1);}
	}
	//-----swap
	public static <T> void swap(T[] a,int l,int r){
		if(l>=0&&l<a.length&&r>=0&&r<a.length){
		T tmp=a[l];
		a[l]=a[r];
		a[r]=tmp;
		}
	}
	public static<T extends Comparable<? super T>> void bbSort(T[] a){
		int o,i;
		for(o=a.length-1;o>1;o--)for(i=0;i<o;i++)
			if(a[i].compareTo(a[i+1])==1){
				swap(a,i,i+1);
			}
	}
	public static<T extends Comparable<? super T>> void selectSort(T[] a){
		int o,i,min;
		for(o=0;o<a.length-1;o++){
			min=o;
			for(i=o+1;i<a.length;i++)
				if(a[i].compareTo(a[i+1])==-1)min=i;
			swap(a,o,min);
		}
	}
	public static <T extends Comparable<? super T>> void insertSort(T[] a){insertSort(a,0,a.length);}
	public static <T extends Comparable<? super T>> void insertSort(T[] a,int l,int r){
		int i,o;
		for(o=l;o<r;o++){
			T tmp=a[o];
			i=o;
			while(i>0&&a[i-1].compareTo(tmp)!=-1){
				a[i]=a[i-1];--i;
			}
			a[i]=tmp;
		}
	}
	public static <T extends Comparable<? super T>> void shellSort(T [] a){
		int i,o;T tmp;
		int h=1;
		while(h<=a.length/3)h=h*3+1;
		while(h>0){
			for(o=h;o<a.length;o++){
				tmp=a[o];
				i=o;
				while(i>h-1&&a[i-h].compareTo(tmp)!=-1){
					a[i]=a[i-h];i-=h;
				}
				a[i]=tmp;
			}//endfor
			h=(h-1)/3;
		}//endwhile
	}
	
	//------------qsort
	public static <T extends Comparable<? super T>> void qsort(T[] a){
		qsort(0,a.length,a);
	}
	private static <T extends Comparable<? super T>> void qsort(int l, int r, T[] a) {
		int size =r-l+1;
		if(size<10) insertSort(a,l,r);
		else{
			int p=partIt(l,r,medianOf3(l,r,a),a);
			qsort(l,p-1,a);
			qsort(p+1,r,a);
		}
	}
	public static <T extends Comparable<? super T>> int partIt(int l,int r,T pivot,T[] a){
		int lp=l,rp=r-1;
		while(true){
			while(a[++lp].compareTo(pivot)==-1);
			while(a[--rp].compareTo(pivot)==1);
			if(lp>=rp)break;
			else swap(a,lp,rp);
		}//endwhile
		swap(a,lp,rp-1);
		return lp;
	}
	private static <T extends Comparable<? super T>> T medianOf3(int l, int r, T[] a) {
		int c=(l+r)/2;
		if(a[l].compareTo(a[c])==1) swap(a,l,c);
		if(a[l].compareTo(a[r])==1) swap(a,l,r);
		if(a[c].compareTo(a[r])==1) swap(a,c,r);
		swap(a,c,r-1);
		return a[r-1];
	}
	//---------mergesort
	public static <T extends Comparable<? super T>>
	void mgsort(T[] a,boolean itr){
		T[] tmp=a.clone();
		if(itr)	mgsort(0,a.length-1,a,tmp);
		else mgsort(a,tmp);
	}
	//-----diedai Merge
	private static <T extends Comparable<? super T>>
	void mgsort(T[] a,T[] tmp){
		int s=1;
		while(s<a.length){
			mgPass(a,tmp,s,a.length);s<<=1;
			mgPass(tmp,a,s,a.length);s<<=1;
		}
	}
	private static <T extends Comparable<? super T>>
	void mgPass(T[] x,T[] y,int s,int n){
		int i=0;
		while(i<=n-(s<<1)){
			mg(i,i+s,i+(s<<1)-1,x,y);
			i+=s<<1;
		}
		if(i+s<n)mg(i,i+s,n-1,x,y);
		else for(int j=i;j<n;j++)y[j]=x[j];
	}
	//------mg
	public static <T extends Comparable<? super T>>
	void mg(int lp,int rp,int re,T[] a,T[] tmp){
		int le=rp-1;int tp=lp,st=lp;
		while(lp<=le&&rp<=re)
			if(a[lp].compareTo(a[rp])<=0)
				tmp[tp++]=a[lp++];
			else
				tmp[tp++]=a[rp++];
		while(lp<=le) tmp[tp++]=a[lp++];
		while(rp<=re) tmp[tp++]=a[rp++];
		for(;re>=st;re--)
			a[re]=tmp[re];
	}
	//-------diguiMerge
	private static <T extends Comparable<? super T>>
	void mgsort(int l,int r,T[] a,T[] tmp){
		if(l<r){
			int c=(l+r)>>>1;
			mgsort(l,c,a,tmp);
			mgsort(c+1,r,a,tmp);
			mg(l,c+1,r,a,tmp);
		}
	}
	public static void swap(int[] a, int l, int r) {
		if(l>=0&&l<a.length&&r>=0&&r<a.length){
			int tmp=a[l];
			a[l]=a[r];
			a[r]=tmp;
			}
		
	}	
}
