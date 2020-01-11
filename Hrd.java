package ran;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Random;
import java.util.TreeMap;


public final class Hrd {
	/**
	 * 2113<br>
	 * 2113<br>
	 * 4665<br>
	 * 4785<br>
	 * 900a
	 */
	
	public static void test(){
		Hrd.setst(null);
		String[] s={"←","→","↑","↓"};
		String[] s1={"","曹操","赵云","黄忠","马超","张飞","关羽","卒甲","卒乙","兵一","兵二"};
		Hrd.print(Hrd.sta);
		long tm =System.currentTimeMillis();
		Hrd.a();
		System.out.println((System.currentTimeMillis()-tm)+":"+Hrd.ans.size());
		for(int i=Hrd.ans.size()-1;i>=0;i--){
			System.out.print(s1[Hrd.ans.get(i)[0]]+s[Hrd.ans.get(i)[1]]+(i%5==4?"\n":","));
		}
	}
	
	
	static String hdlm="2113211346654785900a";
	public static Integer[][] sta=new Integer[7][6];
	public static ArrayList<Integer[]> ans=new ArrayList<Integer[]>();
	public static void a(){
		Queue<Integer[][]> q = new LinkedList<Integer[][]>();
		TreeMap<Integer, Integer[]> m = new TreeMap<Integer, Integer[]>();
		q.add(G.copy(sta));
		m.put(zobhs(sta)[0], null);
		while (!q.isEmpty()) {
			Integer[][] old = q.poll();
			int h1 = zobhs(old)[0];
			int[] h=new int[2];
			if (isWin(old)) {
				ans.clear();
				h[0] = h1;
				while (m.get(h[0]) != null) {
					ans.add(new Integer[] { m.get(h[0])[1], m.get(h[0])[2] });
					h[0] = m.get(h[0])[0];
				}
				print(old);
				return;
			}
			for (int i = 1; i <11; i++)
				for (int j = 0; j<4; j++) {
					sta = G.copy(old);
					if (Move(i, j)) {
						if (isWin(sta)) {
							ans.clear();
							ans.add(new Integer[]{i,j});
							h[0] = h1;
							while (m.get(h[0]) != null) {
								ans.add(new Integer[] { m.get(h[0])[1], m.get(h[0])[2] });
								h[0] = m.get(h[0])[0];
							}
							print(sta);
							return;
						}
						h = zobhs(sta);
						if (!m.containsKey(h[0])) {
							q.add(G.copy(sta));
							m.put(h[0], new Integer[] { h1, i, j });
							m.put(h[1], new Integer[] { h1, i, j });
						}
					}

				}
		}
	}
	
	public static boolean check(){
		int c=1;
		for(int i=1;i<6;i++)for(int j=1;j<5;j++){
			if(sta[i][j]==0)c++;
		}
		System.out.print(c);
		if(c!=2)return false;
		return true;
	}
	
	public static void setst(String s){
		if(s==null)s=hdlm;
		for(int i=0;i<7;i++){
			if(i<6){sta[0][i]=sta[6][i]=0xf;}
			sta[i][0]=sta[i][5]=0xf;
		}
		char[] a=s.toCharArray();
		for(int i=1;i<6;i++)for(int j=1;j<5;j++){
			sta[i][j]=a[(i-1)*4+j-1]=='a'?10:a[(i-1)*4+j-1]-48;
		}	
	}
	
	public static void print(Integer[][] s){
		for(int i=0;i<s.length;i++)for(int j=0;j<s[0].length;j++){
			System.out.print(String.format("%x", s[i][j])+(j%6==5?"\n":","));
		}
	}
	public static int[] type={4,0,2,2,2,2,3,1,1,1,1};
	static int[][][] zob=new int[5][4][5];
	
	static{
		Random rnd=new Random();
		for(int i=0;i<5;i++)for(int j=0;j<4;j++)for(int k=0;k<5;k++)
			zob[i][j][k]=k<4?rnd.nextInt():0;
	}
	public static int zobup(int h,Integer[][] os){
		for(int i=1;i<6;i++)for(int j=1;j<5;j++){
			if(os[i][j]!=sta[i][j]){
				h^=zob[i-1][j-1][type[os[i][j]]];
				h^=zob[i-1][j-1][type[sta[i][j]]];
			}
		}
		return h;
	}
	public static int[] zobhs(Integer[][] sta){
		int h[]=new int[2];
		for(int i=1;i<6;i++)for(int j=1;j<5;j++){
			h[0]^=zob[i-1][j-1][type[sta[i][j]]];
			h[1]^=zob[i-1][4-j][type[sta[i][j]]];
		}
		return h;
	}	
	static int[][] D={{0,-1},{0,1},{-1,0},{1,0}};
	static int[][] type2={{2,2},{1,1},{2,1},{1,2}};
	public static boolean Move(int h,int dir){
		for(int i=1;i<6;i++)for(int j=1;j<5;j++){
			if(sta[i][j]==h&&(sta[i+D[dir][0]][j+D[dir][1]]!=0&&sta[i+D[dir][0]][j+D[dir][1]]!=h))
				return false;
		}
		for(int i=1;i<6;i++)for(int j=1;j<5;j++){
			if(sta[i][j]==h&&sta[i+D[dir][0]][j+D[dir][1]]==0){
				sta[i+D[dir][0]][j+D[dir][1]]=h;
				sta[i+D[dir][0]-type2[type[h]][0]*D[dir][0]][j+D[dir][1]-type2[type[h]][1]*D[dir][1]]=0;
			}
		}
		return true;
	}
	
	public static boolean isWin(Integer[][] sta2){
		return sta2[5][2]==1&&sta2[5][3]==1;
	}
	
	private Hrd(){}
}
