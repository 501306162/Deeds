
/*
	:speech_balloon:		���ƾ������� min-convolution technique...REF:	�� efficient belief propagation for early vision ��
			intoduced in [12]globally optimal registration on a MST using Deeds

	:ballot_box_with_check: 		dt1sq(r+i+k*rl*rl,   indr+i+k*rl*rl,    rl,		-dx,				rl,			    v,z,f, i1);  zoy
			dt1sq(r+j*rl+k*rl*rl,indr+j*rl+k*rl*rl, rl,		-dy,				1,				v,z,f, i1);  zox
			dt1sq(r+i+j*rl,      indr+i+j*rl,       rl,		-dz,				rl*rl,			v,z,f, i1);	xoy
				
			// xxƽ����costֵ   xxƽ��indr�ŵ�  2hw+1  -(������-������)/quant	����ȷ��ƽ���еĵ�
*/
void dt1sq(float *val,             int* ind,     int len,  float offset,        int k,		int* v,float* z,float* f,   int* ind1){
	//
	float INF=1e10;

	int j=0;
	z[0]=-INF;
	z[1]=INF;
	v[0]=0;
	for(int q=1;q<len;q++){
		//	val[q*k]:		ȡ��xxƽ���yy���costֵ
		float s=((val[q*k]+pow((float)q+offset,2.0))-(val[v[j]*k]+pow((float)v[j]+offset,2.0)) ) /(2.0*(float)(q-v[j]));
		while(s<=z[j]){
			j--;
			s=((val[q*k]+pow((float)q+offset,2.0))-(val[v[j]*k]+pow((float)v[j]+offset,2.0)))/(2.0*(float)(q-v[j]));
		}
		j++;
		v[j]=q;
		z[j]=s;
		z[j+1]=INF;
	}
	j=0;
	for(int q=0;q<len;q++){
		f[q]=val[q*k]; //needs to be added to fastDT2 otherwise incorrect
		ind1[q]=ind[q*k];
	} 
	for(int q=0;q<len;q++){
		while(z[j+1]<q){
			j++;
		}
		ind[q*k]=ind1[v[j]];//ind[v[j]*k];
		val[q*k]=pow((float)q-((float)v[j]+offset),2.0)+f[v[j]];//val[v[j]*k];
	}

}


void dt4x(float* r,int* indr,int rl,int lenint,float dx,float dy,float dz,float dq){
	//rl is length of one side
	for(int i=0;i<rl*rl*rl*lenint;i++){
		indr[i]=i;
	}
	int* v=new int[rl]; //slightly faster if not intitialised in each loop
	float* z=new float[rl+1];
	float* f=new float[rl];
	int* i1=new int[rl];
	
	for(int q=0;q<lenint;q++){
		for(int k=0;k<rl;k++){
			for(int i=0;i<rl;i++){
				dt1sq(r+i+k*rl*rl+q*rl*rl*rl,indr+i+k*rl*rl+q*rl*rl*rl,rl,-dx,rl,v,z,f,i1);//);
			}
		}
	}
	for(int q=0;q<lenint;q++){
		for(int k=0;k<rl;k++){
			for(int j=0;j<rl;j++){
				dt1sq(r+j*rl+k*rl*rl+q*rl*rl*rl,indr+j*rl+k*rl*rl+q*rl*rl*rl,rl,-dy,1,v,z,f,i1);//);
			}
		}
	}
	
	for(int q=0;q<lenint;q++){
		for(int j=0;j<rl;j++){
			for(int i=0;i<rl;i++){
				dt1sq(r+i+j*rl+q*rl*rl*rl,indr+i+j*rl+q*rl*rl*rl,rl,-dz,rl*rl,v,z,f,i1);//);
			}
		}
	}

	i1=new int[lenint];
	f=new float[lenint];

	v=new int[lenint];
	z=new float[lenint+1];
	
	for(int k=0;k<rl;k++){
		for(int j=0;j<rl;j++){
			for(int i=0;i<rl;i++){
				dt1sq(r+i+j*rl+k*rl*rl,indr+i+j*rl+k*rl*rl,lenint,-dq,rl*rl*rl,v,z,f,i1);//);
			}
		}
	}
	delete []i1;
	delete []f;

	delete []v;
	delete []z;
	
	i1 = NULL;
	f = NULL;
	v = NULL;
	z = NULL;

	
	
}

/*
		:speech_balloon:	 	dt3x(cost1, inds, len, dx1, dy1, dz1);
				cost1:				 �����ռ���similarity cost
				inds:				 �����ռ��Ǿ���(���)
				len:				 2*hw-1
				dx1,dy1,dz1:		(���ڵ�����-�ӽڵ�����)/quant

		:white_check_mark:		����  fast distance transform

*/

void dt3x(float* r,int* indr,int rl,float dx,float dy,float dz){
	//rl is length of one side
	
	// �����ռ�ĵ�ı���
	for(int i=0;i<rl*rl*rl;i++){
		indr[i]=i;
	}
	//slightly faster if not intitialised in each loop
	int* v=new int[rl]; 
	float* z=new float[rl+1];
	float* f=new float[rl];
	int* i1=new int[rl];
	
	// ��ĳһ��ά���Ͻ��о�������
	for(int k=0;k<rl;k++){
		for(int i=0;i<rl;i++){
			dt1sq(r+i+k*rl*rl,indr+i+k*rl*rl,rl,-dx,rl,v,z,f,i1);// kƬi�� <-> zoy��  <->x��);
		}
	}
	for(int k=0;k<rl;k++){
		for(int j=0;j<rl;j++){
			dt1sq(r+j*rl+k*rl*rl,indr+j*rl+k*rl*rl,rl,-dy,1,v,z,f,i1);//kƬj�� <-> zox��  <->y��);
		}
	}
	
	for(int j=0;j<rl;j++){
		for(int i=0;i<rl;i++){
			dt1sq(r+i+j*rl,indr+i+j*rl,rl,-dz,rl*rl,v,z,f,i1);//i��j�� <-> xoy��  <->z��);
		}
	}
	delete []i1;
	delete []f;

	delete []v;
	delete []z;
	i1 = NULL;
	f = NULL;
	v = NULL;
	z = NULL;
	
}
