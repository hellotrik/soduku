package ran.good;
/*
 * @Prime
 */
public class Prime {
	
	public int p(int n,int[] o){
		int pn=0;
		boolean[] msk=new boolean[n+1];
		for(int i=2;i<n+1;i++)if(!msk[i]){
			o[pn++]=i;
			if(i<32767)for(int j=i*i;j<n+1;j+=i)msk[j]=true;
		}
		return pn;
	}
	private void fp2(long s,int nofyz[],long[] yz,int[] mi){
		if(nofyz[0]==0){
			yz[0]=s;
			mi[0]++;
			nofyz[0]++;
		}else{
			if(yz[nofyz[0]-1]==s)mi[nofyz[0]-1]++;
			else{
					yz[nofyz[0]]=s;
					mi[nofyz[0]++]++;
				}
		}
	}
	public long[] p2(long s){
		if(s<4) return null;
		int[] pr=new int[]{2,3};
		int b=5;
		long[] yz=new long[64];
		int[] mi=new int[64];
		int nofyz[]={0};
		
		int i=0;
		while(s>1){
			for(;i<2;i++)if(s%pr[i]==0){	fp2(pr[i],nofyz,yz,mi);	s/=pr[i];break;	}
			if(i>1){
				for(;b*b<=s;b+=6){
					if(s%b==0){
						fp2(b,nofyz,yz,mi);	s/=b;break;
						}else if(s%(b+2)==0){
							fp2(b+2,nofyz,yz,mi);
							s/=b+2;break;
						}
				}
				if(b*b>s){ fp2(s,nofyz,yz,mi);s=1;}
			}
		}
		
		if(nofyz[0]==1&&mi[0]==1) return null;
		long[] ret=new long[128];
		for(i=0;i<64;i++){
			if(yz[i]>0){
				ret[2*i]=yz[i];
				ret[2*i+1]=mi[i];
			}else break;
		}
		return ret;		
	}
}
