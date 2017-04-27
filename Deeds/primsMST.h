/* Image-driven minimum-spanning-tree calcuation using Prim's algorithm.
 Average run-time should be of n*log(n) complexity.
 Uses heap data structure to speed-up finding the next lowest edge-weight.
 Edge-weights 
*/

class Edge
{
public:
	double weight;
	int vert1;
	int vert2;
	Edge(double w=0,int v1=0,int v2=0);
	bool operator<(const Edge & b) const;
	void print();
};

Edge::Edge(double w,int v1,int v2){
	weight=w;
	vert1=v1;
	vert2=v2;
}
/*
	�� < ����Ϊ >...���ں��������(push_heap,Ĭ��������ͷ��βС)...
					���������Ժ�ͳ���...ͷСβ��
*/
bool Edge::operator<(const Edge & b) const{
	return (this->weight > b.weight);
}


/*
	:speech_balloon:	 last=newEdge(minedge,edgeout,vertices);

		��������ָ���ߵĶ�����յ�...
		�����յ�ֵ,���޸�edgeout...
		�ҵ����Ѿ�����С���ı߼�����ʱ�����յ�ֵΪ-1
*/
int newEdge(Edge edge1,Edge& edgeout,bool* vertices){
	//��С�ߵ� ���������Ƿ񱻷���
	bool new1=vertices[edge1.vert1]; 
	bool new2=vertices[edge1.vert2];
	int out1;
	if(new1^new2){					// 1����2�� ���ҽ���һ�������ʹ�
		if(new1){					// ���1�㱻���ʹ�...����һ������ԭ��(minedge)��ͬ...���Ϊ1,�յ�Ϊ2
			out1=edge1.vert2;
			edgeout=Edge(edge1.weight,edge1.vert1,edge1.vert2);  
		}
		else {						// ���1��û�б����ʹ�...����һ������ԭ�߷����෴...���Ϊ2,�յ�Ϊ1
			out1=edge1.vert1;
			edgeout=Edge(edge1.weight,edge1.vert2,edge1.vert1);
		}
	}
	else{							//�������ʹ�(��û�����ʹ�,Ӧ�ò�����)...���Ѿ�����
		out1=-1;
	}
	return out1;
}

/*
		:speech_balloon:		primsGraph(im1b,ordered1,parents1,step1);
				im1:		ԭʼͼ��(��fixed...Ϊ��)
		:ballot_box_with_check:		ordered:	������ (����ֵ)
				parents:	���ڵ� (����ֵ)
				step:		spacing
*/

void primsGraph(float* im1,int* ordered,int* parents,int step1){
	
	//���ԭʼͼ���ά��С
	int m2=image_m;
	int n2=image_n;
	int o2=image_o;

	//��õ�ǰ������������ά��С
	int m=m2/step1;
	int n=n2/step1;
	int o=o2/step1;
	
	//��ǰ�����������,��Ϊ�ܹ��Ķ�����
	int num_vertices=m*n*o;
	int len=m*n*o;
	timeval time1,time2;

	//�����С  ǰ��..����..����
	int num_neighbours=6;

	//��ǰ��� �ܱ���
	float* edgecost=new float[num_vertices*num_neighbours]; 
	
	//�ܱ�����������
	int* index_neighbours=new int[num_vertices*num_neighbours];
	
	//��ʼ��
	for(int i=0;i<num_vertices*num_neighbours;i++){
		edgecost[i]=0.0;
		index_neighbours[i]=-1;
	}
	
	/*
			���������Ĥmask...
			dx:		2������
			dy:
			dz:
	*/
	int dx[6]={-1,1,0,0,0,0};
	int dy[6]={0,0,-1,1,0,0};
	int dz[6]={0,0,0,0,-1,1};

	int xx,yy,zz,xx2,yy2,zz2;
	
	/*
			calculate edge-weights based on SAD of groups of voxels (for each control-point)
			SAD:	sum of absolute difference
	
	*/
	
	for(int k=0;k<o;k++){		//z
		for(int j=0;j<n;j++){		//x
			for(int i=0;i<m;i++){		//y		�����Ƶ�
				
				
				for(int nb=0;nb<num_neighbours;nb++){		//�������


					/*
							ʹ�����λ��ͼƬ֮��...�жϷ�Χ
							
							[0][0][0]  ֻ�� ��...��...ǰ...(���������)���������
							[.][.][.]  ���������ǵ��Լ���Ե��

							�������ڵĽ��в���...���ڵı���ԭֵ

								D----G
							  /	|	/|		��A��Ϊ��...����a->b Ϊy����
							 A--f C -H						a->c   x����
						   	 |/	  | /						a->d   z����
							 B----E			���������Ű����¹���:
														x����Ϊ0...����Ϊ1
														y	   2...		 3
														z	   4...		 5
					*/
					if((i+dy[nb])>=0&(i+dy[nb])<m&		
						(j+dx[nb])>=0&(j+dx[nb])<n&
						(k+dz[nb])>=0&(k+dz[nb])<o){
						/*	
								���������,û�е����Ϊ-1

								i+j*m+k*m*n		   ��ʾ��ǰ���Ƶ�λ��
								nb*num_vertices    ��ʾ���Ƶ��..���Ϊnb�������
								�洢��ʽΪ:		   һ��num_vertices��С������洢�Ķ������Ϊnb�ŵĶ���...һ��6��
						*/
						index_neighbours[i+j*m+k*m*n+nb*num_vertices]=i+dy[nb]+(j+dx[nb])*m+(k+dz[nb])*m*n;
						//float randv=((float)rand()/float(RAND_MAX));
						//edgecost[i+j*m+k*m*n+nb*num_vertices]=randv;

						/*
							ȷ�� ��ǰ��  ��Ȩֵ

							xx:		ԭʼͼ���п��Ƶ�[xx]λ��...�뵱ǰ���[j]λ�ö�Ӧ
							yy:						[yy]				 [i]	
							zz:						[zz]				 [k]
							
							xx2:	[xx]λ�õ������λ��...
							yy2:	[yy]
							zz2:	[zz]
						*/
						for(int k1=0;k1<step1;k1++){		//z
							for(int j1=0;j1<step1;j1++){		//x
								for(int i1=0;i1<step1;i1++){		//y
									xx=j*step1+j1;			
									yy=i*step1+i1;
									zz=k*step1+k1;

									xx2=(j+dx[nb])*step1+j1;
									yy2=(i+dy[nb])*step1+i1;
									zz2=(k+dz[nb])*step1+k1;
									/*
										 ԭͼ���� ���Ƶ��������� SAD 
										 ��Ϊ ��ǰ��ͼ��� [i][j][k]���Ƶ���[nb]����� ֮��ıߵ�Ȩֵ
									*/
									edgecost[i+j*m+k*m*n+nb*num_vertices]+=fabs(im1[yy+xx*m2+zz*m2*n2]-im1[yy2+xx2*m2+zz2*m2*n2]);
								}
							}
						}
					}
				}
			}
		}
	}
	
	//��ǰ�� ͼ������
	float centrex=n/2;
	float centrey=m/2;
	float centrez=o/2;
	
	//��ǰ�� ͼ������ ���
	int root=m/2+n/2*m+o/2*m*n;
	
	vector<Edge> priority;
	
	// ���п��Ƶ���ʾ���...����������ʱΪtrue,
	bool* vertices=new bool[num_vertices];

	// ���Ƶ� �������ṹ�еĲ㼶����...
	int* level=new int[num_vertices];

	//��ǰ��� ���ڵ�ͷ��ʾ��� ��ʼ��
	for(int i=0;i<num_vertices;i++){
		vertices[i]=false;
		parents[i]=-1;
	}

	//int root=0;

	//���� ���㼶 ��ʼ��
	level[root]=0;

	int last=root;
	//���ڵ������ڵ�����
	vertices[root]=true;
	
	/* 
		Edge(weight,v1,v2)
	*/
	Edge edgeout=Edge(0.0,-1,-1);
	Edge minedge=Edge(0.0,-1,-1);
	float cost=0.0;
	gettimeofday(&time1);
	
	/*
		run n-1 times to have all vertices added
		add edges of new vertex to priority queue	
		������С������,�������Ľڵ���ʾ���,
		ȷ�����Ĳ㼶��ϵ�Լ��Ӹ��ڵ��ϵ,
		�������ıߵ�Ȩֵ�ܺ�,
	*/
	for(int i=0;i<num_vertices-1;i++){ 
		/*
			����ǰ���Ƶ㼰�����б� ����Edge������...ͬʱ����֮ǰ���Ȩֵ
			last+j*num_vertices :  ��ʾ last �ſ��Ƶ�  ��Ӧ�� ��j������
			n :					   ��ʾ last �ſ��Ƶ�  ��Ӧ�� ��j������ ������ͼ����������
		*/
		for(int j=0;j<num_neighbours;j++){
			int n=index_neighbours[last+j*num_vertices];
			//��֤��ǰ���Ƶ��и�����
			if(n>=0){
				priority.push_back(Edge(edgecost[last+j*num_vertices],last,n));
				push_heap(priority.begin(),priority.end());		// �Ѷ�Ϊ��С
			}
		}
		last=-1;
		/*
			find valid edge with lowest weight (main step of Prim's algorithm)
			Prim's algorithm:
						ȷ��һ��������...�����������ĵ㼯��
						������Χ��ѡ��..�����ѡ�߼���...
						�Ӵ�ѡ�߼����������ҵ���СȨֵ�ı�...�����������߼���...
						�����ߵ��յ�..����㼯��...������Ϊ��һ�ε�������
		*/
		while(last==-1){
			minedge=priority.front();				//���ص�һ��Ԫ��...
			pop_heap(priority.begin(),priority.end());
			priority.pop_back();					//ɾ�����һ��Ԫ��
			bool new1=vertices[minedge.vert1];		//is either vertex already part of MST?
			bool new2=vertices[minedge.vert2];
			last=newEdge(minedge,edgeout,vertices); //�޸�edgeout...��Ϊ��minedgeָ�������ı�...�������յ�ֵ
		}
		//��������Ȩֵ��
		cost+=edgeout.weight;
		//���ڵ������С��������
		vertices[last]=true;
		//���ıߵ��յ�ļ��� �����Ļ����ϼ�һ
		level[edgeout.vert2]=level[edgeout.vert1]+1;
		//ȷ�� �����յ�ĸ��ڵ�
		parents[edgeout.vert2]=edgeout.vert1;
	}
	
	/*
	
		find correct ordering in constant time
		
		maxlevel:	��С�����������㼶+1  ===>�ܲ���(����0��)
	
	*/
	int maxlevel=0;
	for(int i=0;i<num_vertices;i++){
		if(level[i]>maxlevel)
			maxlevel=level[i];
	}
	maxlevel++;


	int* leveloffset=new int[maxlevel];
	int* levelcount=new int[maxlevel];
	for(int i=0;i<maxlevel;i++){
		leveloffset[i]=0;
		levelcount[i]=0;
	}
	for(int i=0;i<num_vertices;i++){
		if(level[i]<maxlevel-1)
			
			/*
				counting number of vertices in each level...
				������һ���㼶����(ԭ����0���Ϊ1��)
				����������ڸò�...level[i]+1��...���ж������

				that is, leveloffset ��ʾ���� �ڼ��� �Ķ�����...root��(0��)Ϊ��һ��
			*/
			leveloffset[level[i]+1]++; 

	}
	for(int i=1;i<maxlevel;i++){
		/*
				cumulative sum....
				[1]�㶥����....��Ϊԭ��õ�root�㶥����
				[2]�㶥����....Ϊroot����root+1�㶥�����
				[3]�㶥����....Ϊroot����root+1����root+2�㶥�����
				...
		*/
		leveloffset[i]+=leveloffset[i-1]; 
	}
	/*
		ȷ������˳��ordered
		num:			Ϊ�ܱ���...
		leveloffset:	������ʾ����ƫ����...
						��ǰxxx��Ϊroot�еĵ�...Ȼ����levelcount�����������˳��...��������
	*/
	for(int i=0;i<num_vertices;i++){
		int num=leveloffset[level[i]]+levelcount[level[i]];
		levelcount[level[i]]++;					//counting number of vertices in each level...
		ordered[num]=i;
	}

	gettimeofday(&time2);
	double timeAll=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
	//printf("Prims algorithm with %d levels finished in %f secs.\n",maxlevel,timeAll);
	
	priority.clear();
	
	delete edgecost;
	delete index_neighbours;
	delete levelcount;
	delete leveloffset;
	delete vertices;
	delete level;
	edgecost=NULL;
	index_neighbours = NULL;
	levelcount = NULL;
	leveloffset = NULL;
	vertices = NULL;
	level = NULL;
	
}



