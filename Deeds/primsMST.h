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
	将 < 重载为 >...用于后序堆排序(push_heap,默认排序是头大尾小)...
					重载完了以后就成了...头小尾大
*/
bool Edge::operator<(const Edge & b) const{
	return (this->weight > b.weight);
}


/*
	:speech_balloon:	 last=newEdge(minedge,edgeout,vertices);

		函数用于指定边的顶点和终点...
		返回终点值,并修改edgeout...
		且当边已经在最小树的边集合中时返回终点值为-1
*/
int newEdge(Edge edge1,Edge& edgeout,bool* vertices){
	//最小边的 两个顶点是否被访问
	bool new1=vertices[edge1.vert1]; 
	bool new2=vertices[edge1.vert2];
	int out1;
	if(new1^new2){					// 1点与2点 有且仅有一个被访问过
		if(new1){					// 如果1点被访问过...创建一条边与原边(minedge)相同...起点为1,终点为2
			out1=edge1.vert2;
			edgeout=Edge(edge1.weight,edge1.vert1,edge1.vert2);  
		}
		else {						// 如果1点没有被访问过...创建一条边与原边方向相反...起点为2,终点为1
			out1=edge1.vert1;
			edgeout=Edge(edge1.weight,edge1.vert2,edge1.vert1);
		}
	}
	else{							//都被访问过(都没被访问过,应该不可能)...边已经存在
		out1=-1;
	}
	return out1;
}

/*
		:speech_balloon:		primsGraph(im1b,ordered1,parents1,step1);
				im1:		原始图像(以fixed...为例)
		:ballot_box_with_check:		ordered:	点排序 (返回值)
				parents:	父节点 (返回值)
				step:		spacing
*/

void primsGraph(float* im1,int* ordered,int* parents,int step1){
	
	//获得原始图像各维大小
	int m2=image_m;
	int n2=image_n;
	int o2=image_o;

	//获得当前层的像素数组各维大小
	int m=m2/step1;
	int n=n2/step1;
	int o=o2/step1;
	
	//当前层的总像素数,即为总共的顶点数
	int num_vertices=m*n*o;
	int len=m*n*o;
	timeval time1,time2;

	//邻域大小  前后..左右..上下
	int num_neighbours=6;

	//当前层的 总边数
	float* edgecost=new float[num_vertices*num_neighbours]; 
	
	//总边数排序数组
	int* index_neighbours=new int[num_vertices*num_neighbours];
	
	//初始化
	for(int i=0;i<num_vertices*num_neighbours;i++){
		edgecost[i]=0.0;
		index_neighbours[i]=-1;
	}
	
	/*
			邻域操作掩膜mask...
			dx:		2个邻域
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
			for(int i=0;i<m;i++){		//y		定控制点
				
				
				for(int nb=0;nb<num_neighbours;nb++){		//定邻域点


					/*
							使邻域点位于图片之内...判断范围
							
							[0][0][0]  只有 下...右...前...(按正方向分)三个邻域点
							[.][.][.]  类似其他角点以及边缘点

							在邻域内的进行操作...不在的保持原值

								D----G
							  /	|	/|		以A点为例...坐标a->b 为y正向
							 A--f C -H						a->c   x正向
						   	 |/	  | /						a->d   z正向
							 B----E			邻域点标号序号按如下规则:
														x负向为0...正向为1
														y	   2...		 3
														z	   4...		 5
					*/
					if((i+dy[nb])>=0&(i+dy[nb])<m&		
						(j+dx[nb])>=0&(j+dx[nb])<n&
						(k+dz[nb])>=0&(k+dz[nb])<o){
						/*	
								给邻域点标号,没有点的仍为-1

								i+j*m+k*m*n		   表示当前控制点位置
								nb*num_vertices    表示控制点的..序号为nb的邻域点
								存储方式为:		   一块num_vertices大小的区域存储的都是序号为nb号的顶点...一共6块
						*/
						index_neighbours[i+j*m+k*m*n+nb*num_vertices]=i+dy[nb]+(j+dx[nb])*m+(k+dz[nb])*m*n;
						//float randv=((float)rand()/float(RAND_MAX));
						//edgecost[i+j*m+k*m*n+nb*num_vertices]=randv;

						/*
							确定 当前层  边权值

							xx:		原始图像中控制点[xx]位置...与当前层的[j]位置对应
							yy:						[yy]				 [i]	
							zz:						[zz]				 [k]
							
							xx2:	[xx]位置的邻域点位置...
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
										 原图像中 控制点与邻域点的 SAD 
										 作为 当前层图像的 [i][j][k]控制点与[nb]邻域点 之间的边的权值
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
	
	//当前层 图像中心
	float centrex=n/2;
	float centrey=m/2;
	float centrez=o/2;
	
	//当前层 图像中心 序号
	int root=m/2+n/2*m+o/2*m*n;
	
	vector<Edge> priority;
	
	// 树中控制点访问矩阵...存在于树中时为true,
	bool* vertices=new bool[num_vertices];

	// 控制点 所在树结构中的层级矩阵...
	int* level=new int[num_vertices];

	//当前层的 父节点和访问矩阵 初始化
	for(int i=0;i<num_vertices;i++){
		vertices[i]=false;
		parents[i]=-1;
	}

	//int root=0;

	//树的 根层级 初始化
	level[root]=0;

	int last=root;
	//根节点存在与节点树中
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
		建立最小生成树,更新树的节点访问矩阵,
		确定树的层级关系以及子父节点关系,
		计算树的边的权值总和,
	*/
	for(int i=0;i<num_vertices-1;i++){ 
		/*
			将当前控制点及其所有边 加入Edge队列中...同时给定之前算的权值
			last+j*num_vertices :  表示 last 号控制点  对应的 第j号邻域
			n :					   表示 last 号控制点  对应的 第j号邻域 在整个图像邻域的序号
		*/
		for(int j=0;j<num_neighbours;j++){
			int n=index_neighbours[last+j*num_vertices];
			//保证当前控制点有该邻域
			if(n>=0){
				priority.push_back(Edge(edgecost[last+j*num_vertices],last,n));
				push_heap(priority.begin(),priority.end());		// 堆顶为最小
			}
		}
		last=-1;
		/*
			find valid edge with lowest weight (main step of Prim's algorithm)
			Prim's algorithm:
						确定一个搜索点...放入生成树的点集合
						将其周围待选边..放入待选边集合...
						从待选边集合中搜索找到最小权值的边...放入生成树边集合...
						并将边的终点..放入点集合...将其作为下一次的搜索点
		*/
		while(last==-1){
			minedge=priority.front();				//传回第一个元素...
			pop_heap(priority.begin(),priority.end());
			priority.pop_back();					//删除最后一个元素
			bool new1=vertices[minedge.vert1];		//is either vertex already part of MST?
			bool new2=vertices[minedge.vert2];
			last=newEdge(minedge,edgeout,vertices); //修改edgeout...变为给minedge指定方向后的边...并返回终点值
		}
		//生成树的权值和
		cost+=edgeout.weight;
		//将节点加入最小生成树中
		vertices[last]=true;
		//树的边的终点的级别 在起点的基础上加一
		level[edgeout.vert2]=level[edgeout.vert1]+1;
		//确定 树的终点的父节点
		parents[edgeout.vert2]=edgeout.vert1;
	}
	
	/*
	
		find correct ordering in constant time
		
		maxlevel:	最小生成树的最大层级+1  ===>总层数(包含0级)
	
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
				进行了一个层级提升(原来的0层变为1层)
				如果顶点属于该层...level[i]+1层...进行顶点计数

				that is, leveloffset 表示的是 第几层 的顶点数...root层(0层)为第一层
			*/
			leveloffset[level[i]+1]++; 

	}
	for(int i=1;i<maxlevel;i++){
		/*
				cumulative sum....
				[1]层顶点数....仍为原求得的root层顶点数
				[2]层顶点数....为root层与root+1层顶点个数
				[3]层顶点数....为root层与root+1层与root+2层顶点个数
				...
		*/
		leveloffset[i]+=leveloffset[i-1]; 
	}
	/*
		确定搜索顺序ordered
		num:			为总标序...
		leveloffset:	用来表示顶点偏移量...
						即前xxx块为root中的点...然后用levelcount来控制其访问顺序...依次类推
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



