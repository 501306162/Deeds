/* Dense (stochastic) displacement sampling (deeds)
 for similarity term computation for each node and label.

 Uses a random sub-sampling (defined by global variable RAND_SAMPLES)
 (this step is not described in MICCAI paper)


 Quantization of label space has to be integer or 0.5 (uses trilinear up-sampling)
 =================================================================================
*/


/*
	:speech_balloon:	 	pthread_create( &thread1, NULL, dataCost, (void *) &cosd1);
			���� similarity term

*/
void *dataCost(void *threadarg)
{
	struct cost_data *my_data;
	my_data = (struct cost_data *) threadarg;  //��:(void *) &cosd1....cosd1Ҳ��cost_data���͵�

	// ʹ�ý��ղ����Ժ����ڲ������г�ʼ��
	float* fixed=my_data->im1;
	float* moving=my_data->im1b;
	float* costall=my_data->costall;
	float alpha=my_data->alpha;
	int hw=my_data->hw;
	int step1=my_data->step1;
	float quant=my_data->quant;
	
	// �����ؼ��������Ƿ���  :thinking:
	bool subpixel=false;
	if(quant==0.5)
		subpixel=true;
	
	//  sapcing/(��*search_step_size)
	float alpha1=(float)step1/(alpha*(float)quant);
	
	float randv=((float)rand()/float(RAND_MAX));
	timeval time1,time2;
	// ԭͼ���ά��С
	int m=image_m;
	int n=image_n;
	int o=image_o;
	int sz=m*n*o;
	
	// spacing=step^3
	int step3=pow(step1,3);

	// ��ǰ���ά��С
	int m1=m/step1;
	int n1=n/step1;
	int o1=o/step1;
	int sz1=m1*n1*o1;
	
	
	/* 
		dense displacement space
		len:	 2*r+1
		len4:	(2r+1)^3
		xs:		 ������ռ��xλ��
		ys:
		zs:
		i+j*len+k*len*len:  (2r+1)^3��С����ռ��[i][j][k]�� λ�ƿռ��  
		(j-hw)*quant:   :thinking:���quantΪ0.5...���ʾ  ��ǰ��λ�ü�ȥ�뾶 ��һ��

	*/
	int len=hw*2+1;
	int len4=pow(len,3);
	float* xs=new float[len4];
	float* ys=new float[len4];
	float* zs=new float[len4];
	
	for(int i=0;i<len;i++){  //y
		for(int j=0;j<len;j++){  //x
			for(int k=0;k<len;k++){  //z
					xs[i+j*len+k*len*len]=(float)((j-hw)*quant);
					ys[i+j*len+k*len*len]=(float)((i-hw)*quant);
					zs[i+j*len+k*len*len]=(float)((k-hw)*quant);

			}
		}
	}
	
	int hw2;
	if(subpixel)   //quant==0.5ʱ,subpixel==true
		hw2=hw;		//�뾶����
	else
		hw2=hw*(int)quant;	//�뾶��Ϊԭ���ļ��� ���� ����
	
	//���� �м�movingͼ��
	float* movingi;  

	int mi=m;
	int ni=n;
	int oi=o;
	
	/*
		interpolation with subsampling factor of 2
	*/
	if(subpixel){
		/*
				quant==0.5ʱ,ִ�������ز���
				��ò�ֵ���movingi
				���xs..ys..zs...��ǰ������뾶



				���м�movingͼ��  ��ά����Ϊԭͼ��2����
				x1:		[i][j][k]��...��Ϊ[j]/2...���õ��x�����һ��...
				y1:									�õ��y�����һ��...
				z1:									�õ��z�����һ��...

				x1..y1...z1����ʾ�м�movingͼ�еĵ���ԭʼmoving�е�λ������
		*/
		mi*=2; ni*=2; oi*=2;
		int szi=mi*ni*oi;
		float* x1=new float[szi];
		float* y1=new float[szi];
		float* z1=new float[szi];
		movingi=new float[szi];
		for(int k=0;k<oi;k++){
			for(int j=0;j<ni;j++){
				for(int i=0;i<mi;i++){
					x1[i+j*mi+k*mi*ni]=0.5*(float)j;
					y1[i+j*mi+k*mi*ni]=0.5*(float)i;
					z1[i+j*mi+k*mi*ni]=0.5*(float)k;
				}
			}
		}
		//��movingi�еĵ�ӳ�䵽moving��...��ֵ�������movingi
		interp3(movingi,moving,x1,y1,z1,mi,ni,oi,m,n,o,false);
		delete []x1;
		delete []y1;
		delete []z1;
		x1 = NULL;
		y1 = NULL;
		z1 = NULL;
		/*
			ʹxs�ָ�Ϊ ��ǰx�����ȥ�뾶...
			  ys	
			  zs
		*/
		for(int i=0;i<len4;i++){
			xs[i]*=2.0;           
			ys[i]*=2.0;
			zs[i]*=2.0;
		}
	}
	else{
		/*
			��ִ�������ز���ʱ,
			�м�ͼ��movingi��ԭͼһ��....
				
		*/
		movingi=new float[sz];
		for(int i=0;i<sz;i++){
			movingi[i]=moving[i];
		}
		
	}
		
	int samples=RAND_SAMPLES;  // 64
	bool randommode=samples<pow(step1,3);  // samples < (spacing^3)  ..?
	
	// maxsampΪ  min{samples  ,  spacing������}
	int maxsamp;
	if(randommode){   
		maxsamp=samples;
	}
	else{
		maxsamp=pow(step1,3);
	}

	// �����ռ��
	float* cost1=new float[len4];
	float* costcount=new float[len4];
	// ��ǰ������е㰴12��һ��ֿ�...�ɷֵ�����
	int frac=(int)(sz1/12);

	// spacing/( ��* min{samples  ,  spacing������} )
	float alpha2=alpha1/(float)maxsamp;

	int xx2,yy2,zz2;
	for(int i=0;i<sz1;i++){
		if((i%frac)==0){			// Ϊ0��fracʱ ���
			cout << "x" << flush;	//  cout<<flush��ʾ�������������������ͽ�cout�������������ˢ��
		}
		int z1=i/(m1*n1);           //�õ���...��ǰ������Ƭ
		int x1=(i-z1*m1*n1)/m1;		//			��ǰ��������
		int y1=i-z1*m1*n1-x1*m1;	//			��ǰ��������
		
		z1*=step1;					//�õ���....ԭͼ�е�λ��
		x1*=step1;
		y1*=step1;
		
		bool boundaries=true; //check image boundaries to save min/max computations
		if(subpixel){
			if(x1*2+(step1-1)*2+hw2>=ni|y1*2+(step1-1)*2+hw2>=mi|z1*2+(step1-1)*2+hw2>=oi)
				boundaries=false;
			if(x1*2-hw2<0|y1*2-hw2<0|z1*2-hw2<0)
				boundaries=false;
		}
		
		else{
			if(x1+(step1-1)+hw2>=ni|y1+(step1-1)+hw2>=mi|z1+(step1-1)+hw2>=oi)
				boundaries=false;
			if(x1-hw2<0|y1-hw2<0|z1-hw2<0)
				boundaries=false;
		}
		
		
		for(int l=0;l<len4;l++){
			cost1[l]=0.0;
		}
		
		for(int j1=0;j1<maxsamp;j1++){
			int i1;
			if(randommode)
				//stochastic sampling for speed-up (~8x faster)
				i1=(int)(rand()*pow(step1,3)/float(RAND_MAX));
			else
				i1=j1;
			int zz=i1/(step1*step1);
			int xx=(i1-zz*step1*step1)/step1;
			int yy=i1-zz*step1*step1-xx*step1;
			
			xx+=x1;
			yy+=y1;
			zz+=z1;

			for(int l=0;l<len4;l++){
				if(!(boundaries)){
					if(subpixel){
						xx2=max(min(xx*2+(int)xs[l],ni-1),0);
						yy2=max(min(yy*2+(int)ys[l],mi-1),0);
						zz2=max(min(zz*2+(int)zs[l],oi-1),0);
					}
					else{
						xx2=max(min(xx+(int)(xs[l]),ni-1),0);
						yy2=max(min(yy+(int)(ys[l]),mi-1),0);
						zz2=max(min(zz+(int)(zs[l]),oi-1),0);
					}
				}
				else{
					if(subpixel){
						xx2=xx*2+(int)xs[l];
						yy2=yy*2+(int)ys[l];
						zz2=zz*2+(int)zs[l];
					}
					else{
						xx2=xx+(int)xs[l];
						yy2=yy+(int)ys[l];
						zz2=zz+(int)zs[l];
					}
				}
				//point-wise similarity term (replace if needed, e.g. with pow( ,2.0)
				cost1[l]+=fabs(fixed[yy+xx*m+zz*m*n]-movingi[yy2+xx2*mi+zz2*mi*ni]);

			}
		}
		
		for(int l=0;l<len4;l++){
			costall[i+l*sz1]=alpha2*cost1[l];
		}
	}

	delete []movingi;

	delete []cost1;
	delete []costcount;
	delete []xs;
	delete []ys;
	delete []zs;
	movingi = NULL;
	cost1 = NULL;
	costcount = NULL;
	xs = NULL;
	ys = NULL;
	zs = NULL;
	return NULL;

}

/*
	:white_check_mark: 		warpImage(warped1,im1,ux,vx,wx);   ...
			uxΪͼ��xά�Ĵ�С...
			
			(���������ݸ���)...(����  x��/step  ��õ����ظ���)
*/
template <typename TypeW> //float
void warpImage(TypeW* warped,TypeW* im1,float* u1,float* v1,float* w1){
	//
	int m=image_m; int n=image_n; int o=image_o; int sz=m*n*o;
	//   :warning: warped1=new float[m*n*o];
	interp3(warped,im1,u1,v1,w1,m,n,o,m,n,o,true);
}