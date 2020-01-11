package ran;

/*
 * Copyright (C) 2007 The Android Open Source Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * Matrix math utilities. These methods operate on OpenGL ES format
 * matrices and vectors stored in float arrays.
 * 
 * Matrices are 4 x 4 column-vector matrices stored in column-major
 * order:
 * <pre>
 *  m[offset +  0] m[offset +  4] m[offset +  8] m[offset + 12]
 *  m[offset +  1] m[offset +  5] m[offset +  9] m[offset + 13]
 *  m[offset +  2] m[offset +  6] m[offset + 10] m[offset + 14]
 *  m[offset +  3] m[offset +  7] m[offset + 11] m[offset + 15]
 * </pre>
 * 
 * Vectors are 4 row x 1 column column-vectors stored in order:
 * <pre>
 * v[offset + 0]
 * v[offset + 1]
 * v[offset + 2]
 * v[offset + 3]
 * </pre>
 *
 */
public class Rmath {
	private static float[] mTmpM=new float[16];

	public static void copyV(float[] a,int starta,float[]b,int startb,int len){
		for(int i=0;i<len;i++)a[i+starta]=b[i+startb];
	}
	@SuppressWarnings("unused")public boolean inSquare(float[] s,float[]p){
		float min=s[0],max=s[0];
		for(int i=1;i<4;i++){
			min=Math.min(min,s[i*3]);
			max=Math.max(max,s[i*3]);
		}
		
		float pl=p[0]-max+min;
		int count=0;
		for(int i=0;i<4;i++){
			
		}
		return true;
	}
	
	public static boolean inLine(float[]line,int o,int s,float[]p){//o是起始位 s是步长
		if(!inRect(line,o,s,p))return false;	
		float[][] v=new float[3][3];
		v[0][0]=line[s+o]-line[o];
		v[0][1]=line[1+s+o]-line[1+o];
		v[1][0]=p[0]-line[o];v[1][1]=p[1]-line[1+o];
		Rmath.crossV(v[0],v[1],v[2]);
		if(v[2][2]==0)return true;
		return false;
	}
	public static boolean inRect(float r[],int o,int s,float[]p){
		float xr=Math.signum(p[0]-r[o])*Math.signum(p[0]-r[s+o]);
		float yr=Math.signum(p[1]-r[1+o])*Math.signum(p[1]-r[s+o+1]);
		return xr<=0&&yr<=0;
	}
	public static boolean isRxR(float r1[],int ofs,int strip,float[]r2,int of,int str){
		return Math.max(r1[ofs],r1[ofs+strip])>=Math.min(r2[of],r2[of+str])&&
				Math.max(r2[of],r2[of+str])>=Math.min(r1[ofs],r1[ofs+strip])&&
				Math.max(r1[ofs+1],r1[1+ofs+strip])>=Math.min(r2[1+of],r2[1+of+str])&&
				Math.max(r2[1+of],r2[1+of+str])>=Math.min(r1[1+ofs],r1[1+ofs+strip]);
	}
	public static boolean isLxL(float[]p,int o,int s,float[]q){
		if(!isRxR(p,o,s,q,0,3))return false;
		float[][]pq=new float[6][3];
		pq[0]=Rmath.subV(p,o,q,0,3);//p1-q1
		pq[1]=Rmath.subV(p,o+s,q,0,3);//p2-q1
		pq[2]=Rmath.subV(q,3,q,0,3);//q2-q1
		pq[3]=Rmath.subV(q,0,p,o,3);//q1-p1
		pq[4]=Rmath.subV(q,3,p,o,3);//q2-p1
		pq[5]=Rmath.subV(p,o+s,p,o,3);//p1-p1
		float a=Rmath.crossV(pq[0],pq[2],null);
		float b=Rmath.crossV(pq[1],pq[2],null);
		float c=Rmath.crossV(pq[3],pq[5],null);
		float d=Rmath.crossV(pq[4],pq[5],null);
		return Math.signum(a)*Math.signum(b)<=0&&Math.signum(c)*Math.signum(d)<=0;
	}
	public static boolean inPolygon(float[] py,int s,float[]p){
		int c=0;
		float[] tmp=new float[py.length+s];
		for(int i=0;i<py.length+s;i++){	tmp[i]=py[i%py.length];		}
		float[] ps={p[0],p[1],0,p[0]-2000,p[1],0};
		for(int i=0;i<py.length/s;i++){
			if(Rmath.inLine(tmp,i*s,s,p))return true;
			if(tmp[i*s+1]!=tmp[i*s+s+1]){
				if(tmp[i*s+1]==p[1]&&tmp[i*s+1]>tmp[i*s+s+1])c++;
				else if(tmp[i*s+s+1]==p[1]&&tmp[i*s+s+1]>tmp[i*s+1])c++;
				else{
					if(Rmath.isLxL(tmp,i*s,s,ps))c++;
				}
				
			}
		}
		return (c%2==1);
	}
	public static float[] subV(float[] a,int oa,float[] b,int ob,int s){float[] c=new float[s];	for(int i=0;i<s;i++){c[i]=a[i+oa]-b[i+ob];
	c[i]=(c[i]<0.1&&c[i]>-0.1)?0:c[i];} return c; }	
	public static float[] subV(float[] a,float[] b){float[] c=new float[a.length];	for(int i=0;i<a.length;i++)c[i]=a[i]-b[i]; return c; }
    public static float[] addV(float[] a,float[] b){float[] c=new float[a.length]; 	for(int i=0;i<a.length;i++)c[i]=a[i]+b[i]; return c; }
    public static void scalM(float[] a,float b){for(int i=0;i<a.length;i++)	a[i]*=b;}
    public static void normV(float[] a){
    	double m=1.0/Rmath.length(a[0],a[1],a[2]);
    	for(int i=0;i<a.length;i++){
    		a[i]*=m;
    	}
    }
    public static float crossV(final float[] a, final float[] b, float[] n) {
    	if(n==null)return a[0] * b[1] - a[1] * b[0];
		n[0] = a[1] * b[2] - a[2] * b[1];
		n[1] = a[2] * b[0] - a[0] * b[2];
		n[2] = a[0] * b[1] - a[1] * b[0];
		return n[2];
	}

	public static double dotV(final float[] a, final float[] b) {	
		float res = 0;
		for (int i=0;i<a.length;i++)
		    res += a[i]*b[i];
		return res;
	}
	public static int sign(double x) {	return (x > 0 ? 1 : 0) - (x < 0 ? 1 : 0);	}
	
	public static boolean getRotMatrix(float[] R, float[] I, float[] gravity,
			float[] geomagnetic, float[] ang) {
		float Ax = gravity[0];
		float Ay = gravity[1];
		float Az = gravity[2];
		final float Ex = geomagnetic[0];
		final float Ey = geomagnetic[1];
		final float Ez = geomagnetic[2];
		float Hx = Ey * Az - Ez * Ay;
		float Hy = Ez * Ax - Ex * Az;
		float Hz = Ex * Ay - Ey * Ax;
		final float normH = (float) Math.sqrt(Hx * Hx + Hy * Hy + Hz * Hz);
		if (normH < 0.1f) {
			// device is close to free fall (or in space?), or close to
			// magnetic north pole. Typical values are > 100.
			return false;
		}
		final float invH = 1.0f / normH;
		Hx *= invH;
		Hy *= invH;
		Hz *= invH;
		final float invA = 1.0f / (float) Math
				.sqrt(Ax * Ax + Ay * Ay + Az * Az);
		Ax *= invA;
		Ay *= invA;
		Az *= invA;
		final float Mx = Ay * Hz - Az * Hy;
		final float My = Az * Hx - Ax * Hz;
		final float Mz = Ax * Hy - Ay * Hx;
		if (R != null) {
			if (R.length == 9) {
				R[0] = Hx;	R[1] = Hy;	R[2] = Hz;
				R[3] = Mx;	R[4] = My;	R[5] = Mz;
				R[6] = Ax;	R[7] = Ay;	R[8] = Az;
				if (ang != null) {
					ang[0] = (float) Math.atan2(R[1], R[4]);
					ang[1] = (float) Math.asin(-R[7]);
					ang[2] = (float) Math.atan2(-R[6], R[8]);
				}
			} else if (R.length == 16) {
				R[0] = Hx;	R[1] = Hy;	R[2] = Hz;	R[3] = 0;
				R[4] = Mx;	R[5] = My;	R[6] = Mz;	R[7] = 0;
				R[8] = Ax;	R[9] = Ay;	R[10] = Az;	R[11] = 0;
				R[12] = 0;	R[13] = 0;	R[14] = 0;	R[15] = 1;
				if (ang != null) {
					ang[0] = (float) Math.atan2(R[1], R[5]);
					ang[1] = (float) Math.asin(-R[9]);
					ang[2] = (float) Math.atan2(-R[8], R[10]);
				}
			}
		}
		if (I != null) {
			// compute the inclination matrix by projecting the geomagnetic
			// vector onto the Z (gravity) and X (horizontal component
			// of geomagnetic vector) axes.
			final float invE = 1.0f / (float) Math.sqrt(Ex * Ex + Ey * Ey + Ez
					* Ez);
			final float c = (Ex * Mx + Ey * My + Ez * Mz) * invE;
			final float s = (Ex * Ax + Ey * Ay + Ez * Az) * invE;
			if (I.length == 9) {
				I[0] = 1;	I[1] = 0;	I[2] = 0;
				I[3] = 0;	I[4] = c;	I[5] = s;
				I[6] = 0;	I[7] = -s;	I[8] = c;
			} else if (I.length == 16) {
				I[0] = 1;	I[1] = 0;	I[2] = 0;
				I[4] = 0;	I[5] = c;	I[6] = s;	
				I[8] = 0;	I[9] = -s;	I[10] = c;
				I[3] = I[7] = I[11] = I[12] = I[13] = I[14] = 0;
				I[15] = 1;
			}
		}

		return true;
	}
	
	public static boolean remapCoordinateSystem(float[] inR, int X, int Y,
            float[] outR)
    {
        if (inR == outR) {
            final float[] temp = mTmpM;
            synchronized(temp) {
                // we don't expect to have a lot of contention
                if (remapCoordinateSystemImpl(inR, X, Y, temp)) {
                    final int size = outR.length;
                    for (int i=0 ; i<size ; i++)
                        outR[i] = temp[i];
                    return true;
                }
            }
        }
        return remapCoordinateSystemImpl(inR, X, Y, outR);
    }

    private static boolean remapCoordinateSystemImpl(float[] inR, int X, int Y,
            float[] outR)
    {
        /*
         * X and Y define a rotation matrix 'r':
         *
         *  (X==1)?((X&0x80)?-1:1):0    (X==2)?((X&0x80)?-1:1):0    (X==3)?((X&0x80)?-1:1):0
         *  (Y==1)?((Y&0x80)?-1:1):0    (Y==2)?((Y&0x80)?-1:1):0    (Y==3)?((X&0x80)?-1:1):0
         *                              r[0] ^ r[1]
         *
         * where the 3rd line is the vector product of the first 2 lines
         *
         */

        final int length = outR.length;
        if (inR.length != length)
            return false;   // invalid parameter
        if ((X & 0x7C)!=0 || (Y & 0x7C)!=0)
            return false;   // invalid parameter
        if (((X & 0x3)==0) || ((Y & 0x3)==0))
            return false;   // no axis specified
        if ((X & 0x3) == (Y & 0x3))
            return false;   // same axis specified

        // Z is "the other" axis, its sign is either +/- sign(X)*sign(Y)
        // this can be calculated by exclusive-or'ing X and Y; except for
        // the sign inversion (+/-) which is calculated below.
        int Z = X ^ Y;

        // extract the axis (remove the sign), offset in the range 0 to 2.
        final int x = (X & 0x3)-1;
        final int y = (Y & 0x3)-1;
        final int z = (Z & 0x3)-1;

        // compute the sign of Z (whether it needs to be inverted)
        final int axis_y = (z+1)%3;
        final int axis_z = (z+2)%3;
        if (((x^axis_y)|(y^axis_z)) != 0)
            Z ^= 0x80;

        final boolean sx = (X>=0x80);
        final boolean sy = (Y>=0x80);
        final boolean sz = (Z>=0x80);

        // Perform R * r, in avoiding actual muls and adds.
        final int rowLength = ((length==16)?4:3);
        for (int j=0 ; j<3 ; j++) {
            final int offset = j*rowLength;
            for (int i=0 ; i<3 ; i++) {
                if (x==i)   outR[offset+i] = sx ? -inR[offset+0] : inR[offset+0];
                if (y==i)   outR[offset+i] = sy ? -inR[offset+1] : inR[offset+1];
                if (z==i)   outR[offset+i] = sz ? -inR[offset+2] : inR[offset+2];
            }
        }
        if (length == 16) {
            outR[3] = outR[7] = outR[11] = outR[12] = outR[13] = outR[14] = 0;
            outR[15] = 1;
        }
        return true;
    }

	
	
	 /**
     * Define a viewing transformation in terms of an eye point, a center of
     * view, and an up vector.
     * 
     * @param a a GL10 interface
     * @param eyeX eye point X
     * @param eyeY eye point Y
     * @param eyeZ eye point Z
     * @param centerX center of view X
     * @param centerY center of view Y
     * @param centerZ center of view Z
     * @param upX up vector X
     * @param upY up vector Y
     * @param upZ up vector Z
     */
	public static void gluLookAtM(float[]m,int offset, float eyeX, float eyeY, float eyeZ,
            float centerX, float centerY, float centerZ, float upX, float upY,
            float upZ) {

        // See the OpenGL GLUT documentation for gluLookAt for a description
        // of the algorithm. We implement it in a straightforward way:
        float fx = centerX - eyeX;
        float fy = centerY - eyeY;
        float fz = centerZ - eyeZ;

        // Normalize f
        float rlf = 1.0f / Rmath.length(fx, fy, fz);
        fx *= rlf;   fy *= rlf;   fz *= rlf;
        // compute s = f x up (x means "cross product")
        float sx = fy * upZ - fz * upY;
        float sy = fz * upX - fx * upZ;
        float sz = fx * upY - fy * upX;
        
        // and normalize s
        float rls = 1.0f / Rmath.length(sx, sy, sz);
        sx *= rls;     sy *= rls;      sz *= rls;
        // compute u = s x f
        float ux = sy * fz - sz * fy;
        float uy = sz * fx - sx * fz;
        float uz = sx * fy - sy * fx;

        m[0+offset] = sx;   m[1+offset] = ux;    m[2+offset] = -fx;	
        m[4+offset] = sy;   m[5+offset] = uy;    m[6+offset] = -fy;	
        m[8+offset] = sz;   m[9+offset] = uz;    m[10+offset] = -fz; 
        m[3+offset]=m[7+offset]=m[11+offset]=m[12+offset] = m[13+offset] = m[14+offset] = 0.0f; 
        m[15+offset] = 1.0f;
        translateM(m,offset,-eyeX, -eyeY, -eyeZ);
    }
	 /**
     * Map object coordinates into window coordinates. gluProject transforms the
     * specified object coordinates into window coordinates using model, proj,
     * and view. The result is stored in win.
     * <p>
     * Note that you can use the OES_matrix_get extension, if present, to get
     * the current modelView and projection matrices.
     * 
     * @param objX object coordinates X
     * @param objY object coordinates Y
     * @param objZ object coordinates Z
     * @param model the current modelview matrix
     * @param modelOffset the offset into the model array where the modelview
     *        maxtrix data starts.
     * @param project the current projection matrix
     * @param projectOffset the offset into the project array where the project
     *        matrix data starts.
     * @param view the current view, {x, y, width, height}
     * @param viewOffset the offset into the view array where the view vector
     *        data starts.
     * @param win the output vector {winX, winY, winZ}, that returns the
     *        computed window coordinates.
     * @param winOffset the offset into the win array where the win vector data
     *        starts.
     * @return A return value of GL_TRUE indicates success, a return value of
     *         GL_FALSE indicates failure.
     */
    public static boolean gluProject(float objX, float objY, float objZ,
            float[] model, int modelOffset, float[] project, int projectOffset,
            float[] viewf, int viewOffset, float[] win, int winOffset) {
    	int[] view=new int[4];
    	for(int i=0;i<4;i++)view[i]=(int)viewf[i+viewOffset];
        float[] m = new float[16];
        Rmath.multiplyMM(m, 0, project, projectOffset, model, modelOffset);
        float[] v = new float[4];
        v[0] = objX;
        v[1] = objY;
        v[2] = objZ;
        v[3] = 1.0f;

        float[] v2 = new float[4];       Rmath.multiplyMV(v2, 0, m, 0, v, 0);
        
        float w = v2[3];
        if (w == 0.0f) return false;

        float rw = 1.0f / w;

		win[winOffset] = view[0] + view[2]
				* (v2[0] * rw + 1.0f) * 0.5f;
		win[winOffset + 1] = view[1] + view[3]
				* (v2[1] * rw + 1.0f) * 0.5f;
		win[winOffset + 2] = (v2[2] * rw + 1.0f) * 0.5f;

        return true;
    }

    /**
     * Map window coordinates to object coordinates. gluUnProject maps the
     * specified window coordinates into object coordinates using model, proj,
     * and view. The result is stored in obj.
     * <p>
     * Note that you can use the OES_matrix_get extension, if present, to get
     * the current modelView and projection matrices.
     * 
     * @param winX window coordinates X
     * @param winY window coordinates Y
     * @param winZ window coordinates Z
     * @param model the current modelview matrix
     * @param modelOffset the offset into the model array where the modelview
     *        maxtrix data starts.
     * @param project the current projection matrix
     * @param projectOffset the offset into the project array where the project
     *        matrix data starts.
     * @param view the current view, {x, y, width, height}
     * @param viewOffset the offset into the view array where the view vector
     *        data starts.
     * @param obj the output vector {objX, objY, objZ}, that returns the
     *        computed object coordinates.
     * @param objOffset the offset into the obj array where the obj vector data
     *        starts.
     * @return A return value of GL10.GL_TRUE indicates success, a return value
     *         of GL10.GL_FALSE indicates failure.
     */
    
        
    public static boolean gluUnProject(float winX, float winY, float winZ,
            float[] model, int modelOffset, float[] project, int projectOffset,
            int[] view, int viewOffset, float[] obj, int objOffset) {
        float[] pm = new float[16];
        Rmath.multiplyMM(pm, 0, project, projectOffset, model, modelOffset);

        float[] invPM = new float[16];
        if (!Rmath.invertM(invPM, 0, pm, 0)) {
            return false;
        }

        float[] v = new float[4];

		v[0] = 2.0f * (winX - view[viewOffset + 0]) / view[viewOffset + 2]
				- 1.0f;
		v[1] = 2.0f * (winY - view[viewOffset + 1]) / view[viewOffset + 3]
				- 1.0f;
		v[2] = 2.0f * winZ - 1.0f;
		v[3] = 1.0f;

        float[] v2 = new float[4];

        Rmath.multiplyMV(v2, 0, invPM, 0, v, 0);
        
        float w = v2[3];
        if (w == 0.0f) return false;
        float rw = 1.0f / w;
        obj[objOffset] = v2[0]*rw;
        obj[objOffset + 1] = v2[1]*rw;
        obj[objOffset + 2] = v2[2]*rw;
        return true;
    }

    /**
     * Multiply two 4x4 matrices together and store the result in a third 4x4
     * matrix. In matrix notation: result = lhs x rhs. Due to the way
     * matrix multiplication works, the result matrix will have the same
     * effect as first multiplying by the rhs matrix, then multiplying by
     * the lhs matrix. This is the opposite of what you might expect.
     * 
     * The same float array may be passed for result, lhs, and/or rhs. However,
     * the result element values are undefined if the result elements overlap
     * either the lhs or rhs elements.
     * 
     * @param result The float array that holds the result.
     * @param resultOffset The offset into the result array where the result is
     *        stored.
     * @param lhs The float array that holds the left-hand-side matrix.
     * @param lhsOffset The offset into the lhs array where the lhs is stored
     * @param rhs The float array that holds the right-hand-side matrix.
     * @param rhsOffset The offset into the rhs array where the rhs is stored.
     * 
     * @throws IllegalArgumentException if result, lhs, or rhs are null, or if
     * resultOffset + 16 > result.length or lhsOffset + 16 > lhs.length or
     * rhsOffset + 16 > rhs.length.
     */
	
	public static void multiplyMN(int mid,float[] res,int resoffset,
			float[] lhs,int lhsOffset, float[] rhs, int rhsOffset){
		int lhsrow=(lhs.length-lhsOffset)/mid,rhscol=(rhs.length-rhsOffset)/mid;
		if(res.length-resoffset>=lhsrow*rhscol){
			for(int i=0;i<lhsrow;i++)for(int j=0;j<rhscol;j++){
				res[resoffset+i*rhscol+j]=0;
				for(int k=0;k<mid;k++)
				res[resoffset+i*rhscol+j]+=lhs[lhsOffset+i*mid+k]*rhs[rhsOffset+k*rhscol+j];
			}
		}
	}
	
	
	
    public static void multiplyMM(float[] result, int resultOffset,
            float[] lhs, int lhsOffset, float[] rhs, int rhsOffset){
    	android.opengl.Matrix.multiplyMM(result, resultOffset, lhs, lhsOffset, rhs, rhsOffset);
    };

    /**
     * Multiply a 4 element vector by a 4x4 matrix and store the result in a 4
     * element column vector. In matrix notation: result = lhs x rhs
     * 
     * The same float array may be passed for resultVec, lhsMat, and/or rhsVec.
     * However, the resultVec element values are undefined if the resultVec
     * elements overlap either the lhsMat or rhsVec elements.
     * 
     * @param resultVec The float array that holds the result vector.
     * @param resultVecOffset The offset into the result array where the result
     *        vector is stored.
     * @param lhsMat The float array that holds the left-hand-side matrix.
     * @param lhsMatOffset The offset into the lhs array where the lhs is stored
     * @param rhsVec The float array that holds the right-hand-side vector.
     * @param rhsVecOffset The offset into the rhs vector where the rhs vector
     *        is stored.
     *
     * @throws IllegalArgumentException if resultVec, lhsMat,
     * or rhsVec are null, or if resultVecOffset + 4 > resultVec.length
     * or lhsMatOffset + 16 > lhsMat.length or
     * rhsVecOffset + 4 > rhsVec.length.
     */
    public static void multiplyMV(float[] resultVec,
            int resultVecOffset, float[] lhsMat, int lhsMatOffset,
            float[] rhsVec, int rhsVecOffset){
    	android.opengl.Matrix.multiplyMV(resultVec, resultVecOffset, lhsMat, lhsMatOffset, rhsVec, rhsVecOffset);
    };
    
    /**
     * Transposes a 4 x 4 matrix.
     * 
     * @param mTrans the array that holds the output inverted matrix
     * @param mTransOffset an offset into mInv where the inverted matrix is
     *        stored.
     * @param m the input array
     * @param mOffset an offset into m where the matrix is stored.
     */
    public static void transposeM(float[] mTrans, int mTransOffset, float[] m,
            int mOffset) {
        for (int i = 0; i < 4; i++) {
            int mBase = i * 4 + mOffset;
            mTrans[i + mTransOffset] = m[mBase];
            mTrans[i + 4 + mTransOffset] = m[mBase + 1];
            mTrans[i + 8 + mTransOffset] = m[mBase + 2];
            mTrans[i + 12 + mTransOffset] = m[mBase + 3];
        }
    }

    /**
     * Inverts a 4 x 4 matrix.
     * 
     * @param mInv the array that holds the output inverted matrix
     * @param mInvOffset an offset into mInv where the inverted matrix is
     *        stored.
     * @param m the input array
     * @param mOffset an offset into m where the matrix is stored.
     * @return true if the matrix could be inverted, false if it could not.
     */
    public static boolean invertM(float[] mInv, int mInvOffset, float[] m,
            int mOffset) {
        // Invert a 4 x 4 matrix using Cramer's Rule

        // array of transpose source matrix
        float[] src = new float[16];

        // transpose matrix
        transposeM(src, 0, m, mOffset);

        // temp array for pairs
        float[] tmp = new float[12];

        // calculate pairs for first 8 elements (cofactors)
        tmp[0] = src[10] * src[15];
        tmp[1] = src[11] * src[14];
        tmp[2] = src[9] * src[15];
        tmp[3] = src[11] * src[13];
        tmp[4] = src[9] * src[14];
        tmp[5] = src[10] * src[13];
        tmp[6] = src[8] * src[15];
        tmp[7] = src[11] * src[12];
        tmp[8] = src[8] * src[14];
        tmp[9] = src[10] * src[12];
        tmp[10] = src[8] * src[13];
        tmp[11] = src[9] * src[12];

        // Holds the destination matrix while we're building it up.
        float[] dst = new float[16];

        // calculate first 8 elements (cofactors)
        dst[0] = tmp[0] * src[5] + tmp[3] * src[6] + tmp[4] * src[7];
        dst[0] -= tmp[1] * src[5] + tmp[2] * src[6] + tmp[5] * src[7];
        dst[1] = tmp[1] * src[4] + tmp[6] * src[6] + tmp[9] * src[7];
        dst[1] -= tmp[0] * src[4] + tmp[7] * src[6] + tmp[8] * src[7];
        dst[2] = tmp[2] * src[4] + tmp[7] * src[5] + tmp[10] * src[7];
        dst[2] -= tmp[3] * src[4] + tmp[6] * src[5] + tmp[11] * src[7];
        dst[3] = tmp[5] * src[4] + tmp[8] * src[5] + tmp[11] * src[6];
        dst[3] -= tmp[4] * src[4] + tmp[9] * src[5] + tmp[10] * src[6];
        dst[4] = tmp[1] * src[1] + tmp[2] * src[2] + tmp[5] * src[3];
        dst[4] -= tmp[0] * src[1] + tmp[3] * src[2] + tmp[4] * src[3];
        dst[5] = tmp[0] * src[0] + tmp[7] * src[2] + tmp[8] * src[3];
        dst[5] -= tmp[1] * src[0] + tmp[6] * src[2] + tmp[9] * src[3];
        dst[6] = tmp[3] * src[0] + tmp[6] * src[1] + tmp[11] * src[3];
        dst[6] -= tmp[2] * src[0] + tmp[7] * src[1] + tmp[10] * src[3];
        dst[7] = tmp[4] * src[0] + tmp[9] * src[1] + tmp[10] * src[2];
        dst[7] -= tmp[5] * src[0] + tmp[8] * src[1] + tmp[11] * src[2];

        // calculate pairs for second 8 elements (cofactors)
        tmp[0] = src[2] * src[7];
        tmp[1] = src[3] * src[6];
        tmp[2] = src[1] * src[7];
        tmp[3] = src[3] * src[5];
        tmp[4] = src[1] * src[6];
        tmp[5] = src[2] * src[5];
        tmp[6] = src[0] * src[7];
        tmp[7] = src[3] * src[4];
        tmp[8] = src[0] * src[6];
        tmp[9] = src[2] * src[4];
        tmp[10] = src[0] * src[5];
        tmp[11] = src[1] * src[4];

        // calculate second 8 elements (cofactors)
        dst[8] = tmp[0] * src[13] + tmp[3] * src[14] + tmp[4] * src[15];
        dst[8] -= tmp[1] * src[13] + tmp[2] * src[14] + tmp[5] * src[15];
        dst[9] = tmp[1] * src[12] + tmp[6] * src[14] + tmp[9] * src[15];
        dst[9] -= tmp[0] * src[12] + tmp[7] * src[14] + tmp[8] * src[15];
        dst[10] = tmp[2] * src[12] + tmp[7] * src[13] + tmp[10] * src[15];
        dst[10] -= tmp[3] * src[12] + tmp[6] * src[13] + tmp[11] * src[15];
        dst[11] = tmp[5] * src[12] + tmp[8] * src[13] + tmp[11] * src[14];
        dst[11] -= tmp[4] * src[12] + tmp[9] * src[13] + tmp[10] * src[14];
        dst[12] = tmp[2] * src[10] + tmp[5] * src[11] + tmp[1] * src[9];
        dst[12] -= tmp[4] * src[11] + tmp[0] * src[9] + tmp[3] * src[10];
        dst[13] = tmp[8] * src[11] + tmp[0] * src[8] + tmp[7] * src[10];
        dst[13] -= tmp[6] * src[10] + tmp[9] * src[11] + tmp[1] * src[8];
        dst[14] = tmp[6] * src[9] + tmp[11] * src[11] + tmp[3] * src[8];
        dst[14] -= tmp[10] * src[11] + tmp[2] * src[8] + tmp[7] * src[9];
        dst[15] = tmp[10] * src[10] + tmp[4] * src[8] + tmp[9] * src[9];
        dst[15] -= tmp[8] * src[9] + tmp[11] * src[10] + tmp[5] * src[8];

        // calculate determinant
        float det =
                src[0] * dst[0] + src[1] * dst[1] + src[2] * dst[2] + src[3]
                        * dst[3];

        if (det == 0.0f) {

        }

        // calculate matrix inverse
        det = 1 / det;
        for (int j = 0; j < 16; j++)
            mInv[j + mInvOffset] = dst[j] * det;

        return true;
    }

    /**
     * Computes an orthographic projection matrix.
     * 
     * @param m returns the result
     * @param mOffset
     * @param left
     * @param right
     * @param bottom
     * @param top
     * @param near
     * @param far
     */
    public static void orthoM(float[] m, int mOffset,
        float left, float right, float bottom, float top,
        float near, float far) {
        if (left == right) {
            throw new IllegalArgumentException("left == right");
        }
        if (bottom == top) {
            throw new IllegalArgumentException("bottom == top");
        }
        if (near == far) {
            throw new IllegalArgumentException("near == far");
        }

        final float r_width  = 1.0f / (right - left);
        final float r_height = 1.0f / (top - bottom);
        final float r_depth  = 1.0f / (far - near);
        final float x =  2.0f * (r_width);
        final float y =  2.0f * (r_height);
        final float z = -2.0f * (r_depth);
        final float tx = -(right + left) * r_width;
        final float ty = -(top + bottom) * r_height;
        final float tz = -(far + near) * r_depth;
        m[mOffset + 0] = x;
        m[mOffset + 5] = y;
        m[mOffset +10] = z;
        m[mOffset +12] = tx;
        m[mOffset +13] = ty;
        m[mOffset +14] = tz;
        m[mOffset +15] = 1.0f;
        m[mOffset + 1] = 0.0f;
        m[mOffset + 2] = 0.0f;
        m[mOffset + 3] = 0.0f;
        m[mOffset + 4] = 0.0f;
        m[mOffset + 6] = 0.0f;
        m[mOffset + 7] = 0.0f;
        m[mOffset + 8] = 0.0f;
        m[mOffset + 9] = 0.0f;
        m[mOffset + 11] = 0.0f;
    }
    
    
    /**
     * Define a projection matrix in terms of six clip planes
     * @param m the float array that holds the perspective matrix
     * @param offset the offset into float array m where the perspective
     * matrix data is written
     * @param left
     * @param right
     * @param bottom
     * @param top
     * @param near
     * @param far
     */
    
    public static void frustumM(float[] m, int offset,
            float left, float right, float bottom, float top,
            float near, float far) {
        if (left == right) {
            throw new IllegalArgumentException("left == right");
        }
        if (top == bottom) {
            throw new IllegalArgumentException("top == bottom");
        }
        if (near == far) {
            throw new IllegalArgumentException("near == far");
        }
        if (near <= 0.0f) {
            throw new IllegalArgumentException("near <= 0.0f");
        }
        if (far <= 0.0f) {
            throw new IllegalArgumentException("far <= 0.0f");
        }
        final float r_width  = 1.0f / (right - left);
        final float r_height = 1.0f / (top - bottom);
        final float r_depth  = 1.0f / (near - far);
        final float x = 2.0f * (near * r_width);
        final float y = 2.0f * (near * r_height);
        final float A = 2.0f * ((right + left) * r_width);
        final float B = (top + bottom) * r_height;
        final float C = (far + near) * r_depth;
        final float D = 2.0f * (far * near * r_depth);
        m[offset + 0] = x;
        m[offset + 5] = y;
        m[offset + 8] = A;
        m[offset +  9] = B;
        m[offset + 10] = C;
        m[offset + 14] = D;
        m[offset + 11] = -1.0f;
        m[offset +  1] = 0.0f;
        m[offset +  2] = 0.0f;
        m[offset +  3] = 0.0f;
        m[offset +  4] = 0.0f;
        m[offset +  6] = 0.0f;
        m[offset +  7] = 0.0f;
        m[offset + 12] = 0.0f;
        m[offset + 13] = 0.0f;
        m[offset + 15] = 0.0f;
    }

    /**
     * Computes the length of a vector
     * 
     * @param x x coordinate of a vector
     * @param y y coordinate of a vector
     * @param z z coordinate of a vector
     * @return the length of a vector
     */
    public static float length(float x, float y, float z) {
        return (float) Math.sqrt(x * x + y * y + z * z);
    }
 
    /**
     * Sets matrix m to the identity matrix.
     * @param sm returns the result
     * @param smOffset index into sm where the result matrix starts
     */
    public static void setIdentityM(float[] sm, int smOffset) {
        for (int i=0 ; i<16 ; i++) {
            sm[smOffset + i] = 0;
        }
        for(int i = 0; i < 16; i += 5) {
            sm[smOffset + i] = 1.0f;
        }
    }

    /**
     * Scales matrix  m by sx, sy, and sz, putting the result in sm
     * @param sm returns the result
     * @param smOffset index into sm where the result matrix starts
     * @param m source matrix
     * @param mOffset index into m where the source matrix starts
     * @param x scale factor x
     * @param y scale factor y
     * @param z scale factor z
     */
    public static void scaleM(float[] sm, int smOffset,
            float[] m, int mOffset,
            float x, float y, float z) {
        for (int i=0 ; i<4 ; i++) {
            int smi = smOffset + i;
            int mi = mOffset + i;
            sm[     smi] = m[     mi] * x;
            sm[ 4 + smi] = m[ 4 + mi] * y;
            sm[ 8 + smi] = m[ 8 + mi] * z;
            sm[12 + smi] = m[12 + mi];
        }
    }

    /**
     * Scales matrix m in place by sx, sy, and sz
     * @param m matrix to scale
     * @param mOffset index into m where the matrix starts
     * @param x scale factor x
     * @param y scale factor y
     * @param z scale factor z
     */
    public static void scaleM(float[] m, int mOffset,
            float x, float y, float z) {
        for (int i=0 ; i<4 ; i++) {
            int mi = mOffset + i;
            m[     mi] *= x;
            m[ 4 + mi] *= y;
            m[ 8 + mi] *= z;
        }
    }
    
    /**
     * Translates matrix m by sx, sy, and sz, putting the result in tm
     * @param tm returns the result
     * @param tmOffset index into sm where the result matrix starts
     * @param m source matrix
     * @param mOffset index into m where the source matrix starts
     * @param x translation factor x
     * @param y translation factor y
     * @param z translation factor z
     */
    public static void translateM(float[] tm, int tmOffset,
            float[] m, int mOffset,
            float x, float y, float z) {
        for (int i=0 ; i<4 ; i++) {
            int tmi = tmOffset + i;
            int mi = mOffset + i;
            tm[12 + tmi] = m[mi] * x + m[4 + mi] * y + m[8 + mi] * z +
                m[12 + mi];
        }
    }
 
    /**
     * Translates matrix m by sx, sy, and sz in place.
     * @param m matrix
     * @param mOffset index into m where the matrix starts
     * @param x translation factor x
     * @param y translation factor y
     * @param z translation factor z
     */
    public static void translateM(
            float[] m, int mOffset,
            float x, float y, float z) {
        for (int i=0 ; i<4 ; i++) {
            int mi = mOffset + i;
            m[12 + mi] += m[mi] * x + m[4 + mi] * y + m[8 + mi] * z;
        }
    }
    
    /**
     * Rotates matrix m by angle a (in degrees) around the axis (x, y, z)
     * @param rm returns the result
     * @param rmOffset index into rm where the result matrix starts
     * @param m source matrix
     * @param mOffset index into m where the source matrix starts
     * @param a angle to rotate in degrees
     * @param x scale factor x
     * @param y scale factor y
     * @param z scale factor z
     */
    public static void rotateM(float[] rm, int rmOffset,
            float[] m, int mOffset,
            float a, float x, float y, float z) {
        float[] r = new float[16];
        setRotateM(r, 0, a, x, y, z);
        multiplyMM(rm, rmOffset, m, mOffset, r, 0);
    }
    
    /**
     * Rotates matrix m in place by angle a (in degrees)
     * around the axis (x, y, z)
     * @param m source matrix
     * @param mOffset index into m where the matrix starts
     * @param a angle to rotate in degrees
     * @param x scale factor x
     * @param y scale factor y
     * @param z scale factor z
     */
    public static void rotateM(float[] m, int mOffset,
            float a, float x, float y, float z) {
        float[] temp = new float[32];
        setRotateM(temp, 0, a, x, y, z);
        multiplyMM(temp, 16, m, mOffset, temp, 0);
        System.arraycopy(temp, 16, m, mOffset, 16);
    }
    
    /**
     * Rotates matrix m by angle a (in degrees) around the axis (x, y, z)
     * @param rm returns the result
     * @param rmOffset index into rm where the result matrix starts
     * @param a angle to rotate in degrees
     * @param x scale factor x
     * @param y scale factor y
     * @param z scale factor z
     */
    public static void setRotateM(float[] rm, int rmOffset,
            float a, float x, float y, float z) {
        rm[rmOffset + 3] = 0;
        rm[rmOffset + 7] = 0;
        rm[rmOffset + 11]= 0;
        rm[rmOffset + 12]= 0;
        rm[rmOffset + 13]= 0;
        rm[rmOffset + 14]= 0;
        rm[rmOffset + 15]= 1;
        a *= (float) (Math.PI / 180.0f);
        float s = (float) Math.sin(a);
        float c = (float) Math.cos(a);
        if (1.0f == x && 0.0f == y && 0.0f == z) {
            rm[rmOffset + 5] = c;   rm[rmOffset + 10]= c;
            rm[rmOffset + 6] = s;   rm[rmOffset + 9] = -s;
            rm[rmOffset + 1] = 0;   rm[rmOffset + 2] = 0;
            rm[rmOffset + 4] = 0;   rm[rmOffset + 8] = 0;
            rm[rmOffset + 0] = 1;
        } else if (0.0f == x && 1.0f == y && 0.0f == z) {
            rm[rmOffset + 0] = c;   rm[rmOffset + 10]= c;
            rm[rmOffset + 8] = s;   rm[rmOffset + 2] = -s;
            rm[rmOffset + 1] = 0;   rm[rmOffset + 4] = 0;
            rm[rmOffset + 6] = 0;   rm[rmOffset + 9] = 0;
            rm[rmOffset + 5] = 1;
        } else if (0.0f == x && 0.0f == y && 1.0f == z) {
            rm[rmOffset + 0] = c;   rm[rmOffset + 5] = c;
            rm[rmOffset + 1] = s;   rm[rmOffset + 4] = -s;
            rm[rmOffset + 2] = 0;   rm[rmOffset + 6] = 0;
            rm[rmOffset + 8] = 0;   rm[rmOffset + 9] = 0;
            rm[rmOffset + 10]= 1;
        } else {
            float len = length(x, y, z);
            if (1.0f != len) {
                float recipLen = 1.0f / len;
                x *= recipLen;
                y *= recipLen;
                z *= recipLen;
            }
            float nc = 1.0f - c;
            float xy = x * y;
            float yz = y * z;
            float zx = z * x;
            float xs = x * s;
            float ys = y * s;
            float zs = z * s;       
            rm[rmOffset +  0] = x*x*nc +  c;
            rm[rmOffset +  4] =  xy*nc - zs;
            rm[rmOffset +  8] =  zx*nc + ys;
            rm[rmOffset +  1] =  xy*nc + zs;
            rm[rmOffset +  5] = y*y*nc +  c;
            rm[rmOffset +  9] =  yz*nc - xs;
            rm[rmOffset +  2] =  zx*nc - ys;
            rm[rmOffset +  6] =  yz*nc + xs;
            rm[rmOffset + 10] = z*z*nc +  c;
        }
    }
    
    /**
     * Converts Euler angles to a rotation matrix
     * @param rm returns the result
     * @param rmOffset index into rm where the result matrix starts
     * @param x angle of rotation, in degrees
     * @param y angle of rotation, in degrees
     * @param z angle of rotation, in degrees
     */
    public static void setRotateEulerM(float[] rm, int rmOffset,
            float x, float y, float z) {
        x *= (float) (Math.PI / 180.0f);
        y *= (float) (Math.PI / 180.0f);
        z *= (float) (Math.PI / 180.0f);
        float cx = (float) Math.cos(x);
        float sx = (float) Math.sin(x);
        float cy = (float) Math.cos(y);
        float sy = (float) Math.sin(y);
        float cz = (float) Math.cos(z);
        float sz = (float) Math.sin(z);
        float cxsy = cx * sy;
        float sxsy = sx * sy;
        
        rm[rmOffset + 0]  =   cy * cz;
        rm[rmOffset + 1]  =  -cy * sz;
        rm[rmOffset + 2]  =   sy;
        rm[rmOffset + 3]  =  0.0f;

        rm[rmOffset + 4]  =  cxsy * cz + cx * sz;
        rm[rmOffset + 5]  = -cxsy * sz + cx * cz;
        rm[rmOffset + 6]  =  -sx * cy;
        rm[rmOffset + 7]  =  0.0f;

        rm[rmOffset + 8]  = -sxsy * cz + sx * sz;
        rm[rmOffset + 9]  =  sxsy * sz + sx * cz;
        rm[rmOffset + 10] =  cx * cy;
        rm[rmOffset + 11] =  0.0f;

        rm[rmOffset + 12] =  0.0f;
        rm[rmOffset + 13] =  0.0f;
        rm[rmOffset + 14] =  0.0f;
        rm[rmOffset + 15] =  1.0f;
    }
    private Rmath(){}
}
