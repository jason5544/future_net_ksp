/********** future_net.cpp **********/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <sys/timeb.h>
#include <errno.h>
#include <unistd.h>
#include <signal.h>
#include <iostream>
#include<algorithm>
#include <vector>
#include <queue>
#include <set>

#include "route.h"
#include "lib_record.h"
#include"sys/timeb.h"


using namespace std;

#define INLINE  static __inline
#define PRINT   printf

#define MAX_RECORD  100
#define MAX_LINE_LEN 4000

/******************************franz' code************************/

//Maximum number of edges.
#define MAX_EV (66)
//50

//Maximum number of nodes.
#define MAX_V (777)
//600
const int BigNode=500;
//parameters for debuging
bool debug_showinputstr=true;
bool debug_showbfsqueue=false;
bool debug_showinitmatrix=false;
bool debug_showresultpath=true;
bool debug_showtopoinfo=true;
bool debug_write=true;
bool debug_showbfsfor14or15=false;



//node number of initial graph
int Node4InitInputGragh;
//edge number of initial graph
int Edg4InitInputGragh;

//the initial graph's topological structure(matrix)
int InitGraphMatix[MAX_V][MAX_V];
//the initial graph's topological structure's index(matrix)
int Index4InitGraphEdgeMatix[MAX_V][MAX_V];
//judge the node whether useless after cleaning useless node.
int ExistNodeafterscaledown[MAX_V];

int index_node[MAX_V];
int node_index[MAX_V];

int ENode[MAX_EV];
bool isENode[MAX_V];
int ENodeSize;

int source,destination;
bool beabletofindpath;
bool initsetparavalue=false;


//properties of every nodes
struct degee4node
{
    int outdegreeofNE;
    int outdegreeofE;
    //int indegreeofNE;
    // int indegreeofE;
} degree4node[MAX_V];

//the parameters about SCC
vector<int>SCC_G[MAX_V];
vector<int>SCC_rG[MAX_V];
vector<int>SCC_vs;
bool SCC_used[MAX_V];
int ithnode_ithSCC[MAX_V];
int SCCSize;
int ENodesnumofithSCC[MAX_V];
int SCCDFSAllENode;
bool vs_firstsccdfs;
bool isDFStodestination;
bool used4toposort[MAX_V];
int ithscc_ithtopo[MAX_V];

//the parameters about topological sort of DAG
int indegreeofTopoNode[MAX_V];
vector<int>TopoEdge[MAX_V];
vector<int>rTopoEdge[MAX_V];
int ithtopo_ithscc[MAX_V];
int ithscchaveNE_seq[MAX_V];
int seqithscchaveNE_ithtopo[MAX_V];

bool scc_istoanodeFromDestwithouthaveEN[MAX_V];
bool istoanodeFromDestwithouthaveEN[MAX_V];
bool istoanodeFromDest[MAX_V];

bool debug_isruningsearch_route=true;
//debug function
struct timeb starttime;
struct timeb endtime;


/************************************************************ 我的代码 *********************************************/

/********************** 存储路径模板 *********************/
template<class T>
class PathArc
{
public:
    int arcNum;         //路径中的边编号
    int headNode;	    //边的头节点
    PathArc<T> *preArc;  // 指向路径的前一条边
};

template<class T>
class PathList
{
public:
    T pathWeight;
    PathArc<T> *pathTail;
    PathList(void);
    ~PathList(void);
    bool insert(int &arcNum, T &arcWeight, int &headNode);
    bool isLoop(const int &currentNodeNum);
};

template<class T>
PathList<T>::PathList(void)
{

    pathTail = NULL;
    pathWeight = 0;
}

template<class T>
PathList<T>::~PathList(void)
{
}

// 像路径中插入新的节点
template<class T>
bool PathList<T>::insert(int &arcNum, T &arcWeight, int &headNode)
{
    PathArc<T> *newArc = new PathArc<T>;
    //if (NULL == newArc)
    //{
    //	cerr << "error in creating new path PathArc; " << endl;
    //	return 0;
    //}
    newArc->arcNum = arcNum;
    newArc->preArc = pathTail;
    newArc->headNode = headNode;
    pathTail = newArc;
    pathWeight = arcWeight + pathWeight;
    return 1;
}

// 判断路径是否有环
template<class T>
bool	PathList<T>::isLoop(const int &currentNodeNum)
{
    //cout << currentNodeNum << " ";
    //for (PathArc<T> *arc = this->pathTail; arc != NULL; arc = arc->preArc)
    //{
    //	cout << arc->headNode << " ";
    //}
    //cout << endl;
    if (this->pathTail == NULL)
    {
        return 0;
    }
    for (PathArc<T> *arc = this->pathTail; arc != NULL; arc = arc->preArc)
    {
        if (arc->headNode == currentNodeNum)
        {
            return 1;
        }
    }

    return 0;
}
/**************************************************************/



/********************** 存储迭代节点模板 *********************/
template<class T>
class KspIterNode
{
public:
    PathList<T> path;		// 路径
    int nodeNum;            // 节点数
    KspIterNode(void);
    KspIterNode(const T &weight,const int &nodeNum);
    ~KspIterNode(void);
    bool operator<(const KspIterNode<T> &node) const;
};

template<class T>
KspIterNode<T>::KspIterNode(void)
{ }

template<class T>
KspIterNode<T>::KspIterNode(const T &weight,const int &nodeNum)
{
    this->nodeNum = nodeNum;
    path.pathWeight = weight;
}

template<class T>
KspIterNode<T>::~KspIterNode(void)
{ }

// 重载 小于号  ,为了应用优先级队列
template<class T>
bool KspIterNode<T>::operator<(const KspIterNode<T> &node)	const
{
    return  this->path.pathWeight > node.path.pathWeight;
}
/**********************************************************************/


/*************************** 存储图 出边表 入边表 ***********************/
class GraphArc
{
public:
    int weight;                //图的边的权重
    int tailNode;				// 边的尾节点
    int next;					// 指向相同出发节点的 下一条边
};


class InGraphArc
{
public:
    int headNode;
    int next;
    //int weight;
};

class InGraph
{
public:
    int *node;
    InGraphArc *inArc;
};


class Graph
{
public:
    int sourse;						//起始节点
    int destination;				//终点
    vector<int> includingSet;		// 必须经过的节点
//	int nodeNum;                    //图的节点数量
    int edgeNum;					//图的边的数量
    int *node;						//存储图的节点数组
    GraphArc *arc;					//存储图的边的数组

    // 字符串转化为数字
    static void strToNum(int &num, char *str,const int &pos0, const int &pos1);
    // 读入图
    static void readGraph(char *graph[5000], const int &edge_num, Graph &linkGraph, InGraph &inLinkGraph);
    // 读入 condition
    static void readCondition(char *condition, Graph &linkGraph);
    // 读图判定图的节点数目
    static int readNodeNum(char *graph[5000], const int &edge_num);
    // 读必过节点数目
    static int readConditionSize(char *condition);

    // K短路算法
    static void kShortestPath (vector<int> &result, const Graph &linkGraph);
    // 预处理图
    static void preprocess(Graph &linkGraph, InGraph &inLinkGraph);

};

// 字符串转化为数字
void Graph::strToNum(int &num, char *str, const int &pos0, const int &pos1)
{
    char *tmp = new char[4];
    int i = 0;
    int j = pos0 + 1;
    for(; j < pos1; )
        tmp[i++]=str[j++];
    tmp[i] = 0;

    num = atoi(tmp);
//	cout << "num = " << num << endl;
    delete[] tmp;
}

// 读图判定图的节点数目
int Graph::readNodeNum(char *graph[5000], const int &edge_num)
{
    set<int> s;
    int *pos = new int[3];
    int edgeNum = edge_num - 1; // 最大边编号
    for (int i = 0; i < edgeNum; i++)
    {
        int j = 0;
        int tag = 0;
        while('\n' != graph[i][j])
        {
            if (',' == graph[i][j])
            {
                pos[tag++] = j;
            }
            j++;
        }
        pos[tag] = j;
        int head, tail;
        Graph::strToNum(head, graph[i], pos[0], pos[1]);
        Graph::strToNum(tail, graph[i], pos[1], pos[2]);
        s.insert(head);
        s.insert(tail);


    }

    int j = 0;
    int tag = 0;
    while(0 != graph[edgeNum][j])
    {
        if (',' == graph[edgeNum][j])
        {
            pos[tag++] = j;
        }
        j++;
    }
    pos[tag] = j;
    int head, tail ;
    Graph::strToNum(head, graph[edgeNum], pos[0], pos[1]);
    Graph::strToNum(tail, graph[edgeNum], pos[1], pos[2]);
    s.insert(head);
    s.insert(tail);

    return s.size();
    delete[] pos;

}

// 读入图,分别用node存储节点,arc 存储边,linkGraph 表示入边表,inLinkGraph 表示出边表,都是基于头插法建立图
void Graph::readGraph(char *graph[5000], const int &edge_num, Graph &linkGraph, InGraph &inLinkGraph)
{
    linkGraph.edgeNum = edge_num;
    linkGraph.arc = new GraphArc[edge_num];
    inLinkGraph.inArc = new InGraphArc[edge_num];

    linkGraph.node = new int[600];
    inLinkGraph.node = new int[600];
    for (int i = 0; i < 600; i++)
    {
        linkGraph.node[i] = -1;
        inLinkGraph.node[i] = -1;
    }

    int *pos = new int[3];
    int edgeNum = linkGraph.edgeNum - 1; // 最大边编号
    for (int i = 0; i < edgeNum; i++)
    {
        int j = 0;
        int tag = 0;
        while('\n' != graph[i][j])
        {
            if (',' == graph[i][j])
            {
                pos[tag++] = j;
            }
            j++;
        }
        pos[tag] = j;
        int head, tail, weight;
        Graph::strToNum(head, graph[i], pos[0], pos[1]);
        Graph::strToNum(tail, graph[i], pos[1], pos[2]);
        Graph::strToNum(weight, graph[i], pos[2], pos[3]);
        linkGraph.arc[i].weight = weight;
        //inLinkGraph.inArc[i].weight = weight;
        inLinkGraph.inArc[i].headNode = head;
        linkGraph.arc[i].tailNode = tail;

        linkGraph.arc[i].next = linkGraph.node[head];
        inLinkGraph.inArc[i].next = inLinkGraph.node[tail];
        linkGraph.node[head] = i;
        inLinkGraph.node[tail] = i;
    }

    //  topo  最后一行 没有 回车, 所以单独读入
    int j = 0;
    int tag = 0;
    while(0 != graph[edgeNum][j])
    {
        if (',' == graph[edgeNum][j])
        {
            pos[tag++] = j;
        }
        j++;
    }
    pos[tag] = j;
    int head, tail , weight;
    Graph::strToNum(head, graph[edgeNum], pos[0], pos[1]);
    Graph::strToNum(tail, graph[edgeNum], pos[1], pos[2]);
    Graph::strToNum(weight, graph[edgeNum], pos[2], pos[3]);
    inLinkGraph.inArc[edgeNum].headNode = head;
    linkGraph.arc[edgeNum].tailNode = tail;
    linkGraph.arc[edgeNum].weight = weight;
    //inLinkGraph.inArc[edgeNum].weight = weight;
    linkGraph.arc[edgeNum].next = linkGraph.node[head];
    inLinkGraph.inArc[edgeNum].next = inLinkGraph.node[tail];
    linkGraph.node[head] = edgeNum;
    inLinkGraph.node[tail] = edgeNum;

    delete[] pos;
}

// 计算必过节点数目
int Graph::readConditionSize(char *condition)
{
    int i = 0;
    int *pos = new int[2];
    int tag = 0;
    while(tag < 2)
    {
        if (',' == condition[i])
        {
            pos[tag++] = i;
        }
        i++;
    }
    //	cout << "pos[0] = " << pos[0] << " " << " pos[1] = " << pos[1] << endl;
    //Graph::strToNum(linkGraph.sourse, condition, -1, pos[0]);
    //Graph::strToNum(linkGraph.destination, condition, pos[0], pos[1]);
    delete[] pos;

    set<int> s;
    int *pos1 = new int[50];

    int tag1 = 0;
    pos1[tag1++] = i - 1;
    while(0 != condition[i])
    {
        if ('|' == condition[i])
        {
            pos1[tag1++] = i;
        }
        i++;
    }
    pos1[tag1] = i;

    int includeNode;
    for (int i = 0; i < tag1; i++)
    {
        //cout << "pos1[i] = " << pos1[i] << endl;
        Graph::strToNum(includeNode, condition, pos1[i], pos1[i + 1]);
        //cout << "includeNode = " << includeNode << endl;
        s.insert(includeNode);
    }

    return s.size();
    delete[] pos1;
}

// 读入必过节点 condition
void Graph::readCondition(char *condition, Graph &linkGraph)
{
    int i = 0;
    int *pos = new int[2];
    int tag = 0;
    while(tag < 2)
    {

        if (',' == condition[i])
        {
            pos[tag++] = i;
        }
        i++;
    }
    //	cout << "pos[0] = " << pos[0] << " " << " pos[1] = " << pos[1] << endl;
    Graph::strToNum(linkGraph.sourse, condition, -1, pos[0]);
    Graph::strToNum(linkGraph.destination, condition, pos[0], pos[1]);
    delete[] pos;

    int *pos1 = new int[50];

    int tag1 = 0;
    pos1[tag1++] = i - 1;
    while(0 != condition[i])
    {
        if ('|' == condition[i])
        {
            pos1[tag1++] = i;
        }
        i++;
    }
    pos1[tag1] = i;

    int includeNode;
    for (int i = 0; i < tag1; i++)
    {
        //cout << "pos1[i] = " << pos1[i] << endl;
        Graph::strToNum(includeNode, condition, pos1[i], pos1[i + 1]);
        //cout << "includeNode = " << includeNode << endl;
        linkGraph.includingSet.push_back(includeNode);
    }
    delete[] pos1;
}

// 对图做预处理
void Graph::preprocess(Graph & linkGraph, InGraph &inLinkGraph)
{
    int minWeight = 13;
    int i = inLinkGraph.node[linkGraph.destination];

    i = inLinkGraph.node[linkGraph.destination];
    while(i > -1)
    {
        linkGraph.arc[i].weight = linkGraph.arc[i].weight - minWeight;
        i = inLinkGraph.inArc[i].next;
    }

    for(int i = 0; i < linkGraph.includingSet.size(); i++)
    {

        int j = linkGraph.node[linkGraph.includingSet[i]];

        while(j > -1)
        {
            linkGraph.arc[j].weight = linkGraph.arc[j].weight - minWeight;
            j = linkGraph.arc[j].next;
        }

        j = inLinkGraph.node[linkGraph.includingSet[i]];
        while(j > -1)
        {
            linkGraph.arc[j].weight = linkGraph.arc[j].weight - minWeight;
            j = inLinkGraph.inArc[j].next;
        }
    }
}


// K短路算法
void Graph::kShortestPath (vector<int> &result , const Graph &linkGraph)
{
    priority_queue< KspIterNode<int> > ahpHeap;

    KspIterNode<int> currentNode(0, linkGraph.sourse);
    ahpHeap.push(currentNode);

    int k = 0;
    int minCost = 999;
    PathArc<int> *minPath = NULL;
    while (!ahpHeap.empty() && k < 20)
    {
        currentNode = ahpHeap.top();
        ahpHeap.pop();

        if (currentNode.nodeNum == linkGraph.destination)
        {
            int i = 0;
            for(; i < linkGraph.includingSet.size(); i++)
            {
                PathArc<int> *pathArc = currentNode.path.pathTail;
                for(; pathArc != NULL; pathArc = pathArc->preArc)
                {
                    if(pathArc->headNode == linkGraph.includingSet[i])
                    {
                        break;
                    }

                }
                if(NULL == pathArc)
                {
                    break;
                }
            }
            if(i == linkGraph.includingSet.size())
            {
                if(currentNode.path.pathWeight < minCost)
                {
                    minCost = currentNode.path.pathWeight;
                    minPath = currentNode.path.pathTail;
                }
                //	cout << "迭代出" << k << "条路径" << " weight = " << currentNode.path.pathWeight << endl;
                k++;
            }
            else
            {
                delete currentNode.path.pathTail;
            }
            continue;
        }

        for (int i = linkGraph.node[currentNode.nodeNum]; i > -1; i = linkGraph.arc[i].next)
        {
            if (currentNode.path.isLoop(linkGraph.arc[i].tailNode))
            {
                continue;
            }
            KspIterNode<int> nextNode;
            nextNode.nodeNum = linkGraph.arc[i].tailNode;
            nextNode.path = currentNode.path;
            nextNode.path.insert(i, linkGraph.arc[i].weight, currentNode.nodeNum);
            ahpHeap.push(nextNode);
        }
    }

    if(minPath != NULL)
    {
        for(PathArc<int> *pathArc = minPath; pathArc != NULL; pathArc = pathArc->preArc)
        {
            result.push_back(pathArc->arcNum);
        }
    }

}

/**********************************************************************/

/**********************************franz's code**************************/

void DebugPrint()
{
    cout<<"source node:"<<node_index[source]<<"  destination node:"<<node_index[destination]<<endl;
    cout<<"essential node"<<ENodeSize<<":";
    for(int i=0; i<ENodeSize; i++)
    {
        cout<<node_index[ENode[i]]<<" ";
    }
    cout<<endl;
    cout<<"node size:"<<Node4InitInputGragh<<endl;
    cout<<"edge size:"<<Edg4InitInputGragh<<endl;

    if(debug_showinitmatrix)
    {
        cout<<"Initmatrix:"<<endl;
        for(int i=0; i<Node4InitInputGragh; i++)
        {
            cout<<node_index[i]<<" ";
        }
        cout<<endl;
        cout<<endl;

        for(int i=0; i<Node4InitInputGragh; i++)
        {
            for(int j=0; j<Node4InitInputGragh; j++)
            {
                cout<<InitGraphMatix[i][j]<<" ";
            }
            cout<<endl;
        }
    }

    cout<<"SCC     "<<SCCSize<<"   :";
    for(int i=0; i<Node4InitInputGragh; i++)
    {
        cout<<ithnode_ithSCC[i]<<" ";
    }
    cout<<endl;
    cout<<"ithSCChaveNE"<<":";
    for(int i=0; i<Node4InitInputGragh; i++)
    {
        cout<<ENodesnumofithSCC[ithnode_ithSCC[i]]<<" ";
    }
    cout<<endl;
}
//wrap some node whose indgree lower than 1 or outdgree than 1.
void scaledowngragh()
{
    bool judgeindegree,judgeoutdegree;
    judgeindegree=judgeoutdegree=false;
    for(int k=2; k<Node4InitInputGragh; k++)
    {
        for(int i=0; i<Node4InitInputGragh; i++)
        {
            //attention
            if(0<InitGraphMatix[k][i])
                judgeoutdegree=true;
        }
        for(int j=0; j<Node4InitInputGragh; j++)
        {
            if(0<InitGraphMatix[j][k])
                judgeindegree=true;
        }
        if((false==judgeindegree)||(false==judgeoutdegree))
        {
            ExistNodeafterscaledown[k]=0;
        }
    }
}
//judge the graph whether succeed to find path or not.
bool judgegraph()
{
    bool judgeindegree,judgeoutdegree;
    judgeindegree=judgeoutdegree=false;
    for(int j=0; j<Node4InitInputGragh; j++)
    {
        if(0<InitGraphMatix[source][j])
            judgeoutdegree=true;
    }
    if(false==judgeoutdegree)
        return false;
    for(int i=0; i<Node4InitInputGragh; i++)
    {
        if(0<InitGraphMatix[i][destination])
            judgeindegree=true;
    }
    if(false==judgeindegree)
        return false;
    for(int k=0; k<ENodeSize; k++)
    {
        if(0==ExistNodeafterscaledown[ENode[k]])
        {
            return false;
        }
    }
    return true;
}
//before finding path clear and set the parameter's value
void initcleardata()
{
    Node4InitInputGragh=0;
    Edg4InitInputGragh=0;
    memset(index_node,-1,sizeof(index_node));
    memset(node_index,-1,sizeof(node_index));
    if(initsetparavalue)
    {
        memset(InitGraphMatix,0,sizeof(InitGraphMatix));
        memset(Index4InitGraphEdgeMatix,0,sizeof(Index4InitGraphEdgeMatix));
        memset(ExistNodeafterscaledown,0,sizeof(ExistNodeafterscaledown));
        memset(isENode,0,sizeof(isENode));
    }
    return ;

}
//add edge to the graph(vector)
void add_edge(int from,int to)
{
    SCC_G[from].push_back(to);
    SCC_rG[to].push_back(from);
    if(true==isENode[to])
    {
        degree4node[from].outdegreeofNE++;
    }
    else
    {
        degree4node[from].outdegreeofE++;
    }
}

//change the physical structure of the graph(matrix->list)
void convertMatrixtoList()
{
    if(initsetparavalue)
    {
        for(int i=0; i<Node4InitInputGragh; i++)
        {
            if(ExistNodeafterscaledown[i]>0)
            {
                SCC_G[i].clear();
            }
        }
    }
    for(int i=0; i<Node4InitInputGragh; i++)
    {
        if(ExistNodeafterscaledown[i]>0)
        {
            for(int j=0; j<Node4InitInputGragh; j++)
            {
                if(i!=j)
                    if(ExistNodeafterscaledown[j]>0)
                    {
                        if(InitGraphMatix[i][j]>0)
                        {
                            add_edge(i,j);
                        }
                    }
            }
        }
    }
}
//forward dfs function for SCC
void ccdfs(int v)
{
    if(true==vs_firstsccdfs)
    {
        used4toposort[v]=true;
        if(true==isENode[v])
            SCCDFSAllENode++;
        if(v==destination)
        {
            isDFStodestination=true;
        }
    }
    SCC_used[v]=true;
    for(int i=0; i<SCC_G[v].size(); i++)
    {
        if(!SCC_used[SCC_G[v][i]])
            ccdfs(SCC_G[v][i]);
    }
    SCC_vs.push_back(v);
}
//reverse dfs function for SCC
void ccrdfs(int v,int k)
{
    SCC_used[v]=true;
    ithnode_ithSCC[v]=k;
    for(int i=0; i<SCC_rG[v].size(); i++)
    {
        if(!SCC_used[SCC_rG[v][i]]) ccrdfs(SCC_rG[v][i],k);
    }
}
//main function for SCC
int scc()
{
    vs_firstsccdfs=true;
    isDFStodestination=false;
    SCCDFSAllENode=0;
    if(true==initsetparavalue)
    {
        memset(SCC_used,0,sizeof(SCC_used));
        memset(used4toposort,0,sizeof(used4toposort));
        SCC_vs.clear();
    }
    for(int v=0; v<Node4InitInputGragh; v++)
    {
        //if(existNode[v]>0)
        if(!SCC_used[v])ccdfs(v);
        vs_firstsccdfs=false;
    }
    memset(SCC_used,0,sizeof(SCC_used));
    int k=0;
    for(int i=SCC_vs.size()-1; i>=0; i--)
    {
        // if(existNode[vs[i]]>0)
        if(!SCC_used[SCC_vs[i]])ccrdfs(SCC_vs[i],k++);
    }
    if(true==initsetparavalue)
        memset(ENodesnumofithSCC,0,sizeof(ENodesnumofithSCC));
    for(int i=0; i<ENodeSize; i++)
    {
        ENodesnumofithSCC[ithnode_ithSCC[ENode[i]]]++;
    }
    return k;
}
//topo sort function
bool TopoOrder(int topon)
{
    int i,top=-1;
    int sorti=0;
    for(i=(topon-1); i>=0; i--)
    {
        if(0==indegreeofTopoNode[i])
        {
            indegreeofTopoNode[i]=top;
            top=i;
        }
    }
    if(debug_showtopoinfo)
        cout<<"TopoSort "<<topon<<"  :";
    for(i=0; i<topon; ++i)
    {
        if(top==-1)
        {
            if(debug_write)
                cout<<"The graph is not DAG"<<endl;
            return false;
        }
        int j=top;
        ithtopo_ithscc[sorti]=j;
        if(debug_showtopoinfo)
            cout<<j<<" ";
        sorti++;
        top=indegreeofTopoNode[top];
        int topofirstenode=false;
        int exmid=0;
        for(int k=0; k<TopoEdge[j].size(); ++k)
        {
            if((0==(--indegreeofTopoNode[TopoEdge[j][k]])))
            {
                if(false==topofirstenode)
                {
                    indegreeofTopoNode[TopoEdge[j][k]]=top;
                    top=TopoEdge[j][k];
                    topofirstenode=true;
                }
                else
                {
                    if(ENodesnumofithSCC[top]>0)
                    {
                        exmid=indegreeofTopoNode[top];
                        indegreeofTopoNode[TopoEdge[j][k]]=exmid;
                        indegreeofTopoNode[top]=TopoEdge[j][k];
                    }
                    else
                    {
                        indegreeofTopoNode[TopoEdge[j][k]]=top;
                        top=TopoEdge[j][k];
                    }
                }
            }
        }
    }
    if(debug_showtopoinfo)
        cout<<endl;
}
//after SCC the function that get the topo sequence.
void getTopoSort4SCC()
{
    bool vs_matrixtopograph[MAX_V][MAX_V];
    int toponodesize=SCCSize;
    if(true==initsetparavalue)
    {
        memset(TopoEdge,0,sizeof(TopoEdge));
        memset(indegreeofTopoNode,0,sizeof(indegreeofTopoNode));
    }
    for(int i=0; i<Node4InitInputGragh; i++)
    {
        if(true==used4toposort[i])
        {
            for(int j=0; j<Node4InitInputGragh; j++)
            {
                if(i!=j)
                    if(true==used4toposort[j])
                    {
                        if(InitGraphMatix[i][j]>0)
                        {
                            if((ithnode_ithSCC[i]!=ithnode_ithSCC[j])&&(false==vs_matrixtopograph[ithnode_ithSCC[i]][ithnode_ithSCC[j]]))//&&(topoedge[ithnode_ithSCC[i]][ithnode_ithSCC[j]]==false))
                            {
                                TopoEdge[ithnode_ithSCC[i]].push_back(ithnode_ithSCC[j]);//=true;
                                rTopoEdge[ithnode_ithSCC[j]].push_back(ithnode_ithSCC[i]);
                                vs_matrixtopograph[ithnode_ithSCC[i]][ithnode_ithSCC[j]]=true;
                                indegreeofTopoNode[ithnode_ithSCC[j]]++;
                            }
                        }
                    }
            }
        }
    }
    if(debug_showtopoinfo)
    {
        cout<<"topoindgree:";
        for(int i=0; i<toponodesize; i++)
        {
            cout<<indegreeofTopoNode[i]<<" ";
        }
        cout<<endl;
        for(int i=0; i<toponodesize; i++)
        {
            for(int j=0; j<TopoEdge[i].size(); j++)
            {
                cout<<TopoEdge[i][j]<<" ";
            }
            cout<<endl;
        }
    }
    TopoOrder(toponodesize);
    for(int i=0; i<toponodesize; i++)
    {
        ithscc_ithtopo[ithtopo_ithscc[i]]=i;
    }
    seqithscchaveNE_ithtopo[0]=ithscc_ithtopo[source];
    //seqnumofitscchaveNE[ithnode_ithSCC[source]]=0;
    //seqnumofitscchaveNE[ithnode_ithSCC[destination]]=seq;
    int seq=1;
    for(int i=0; i<toponodesize; i++)
    {
        if(ENodesnumofithSCC[ithtopo_ithscc[i]]>0)
        {
            ithscchaveNE_seq[ithtopo_ithscc[i]]=seq;
            seqithscchaveNE_ithtopo[seq]=i;//ithtopo_ithscc[ithnode_ithSCC[i]];
            seq++;
        }
    }
    ithscchaveNE_seq[ithnode_ithSCC[destination]]=seq;
    seqithscchaveNE_ithtopo[seq]=ithscc_ithtopo[ithnode_ithSCC[destination]];
    seq++;
    if(debug_showtopoinfo)
    {
        cout<<"toposort:";
        for(int i=0; i<toponodesize; i++)
        {
            cout<<ithtopo_ithscc[i]<<" ";
        }
        cout<<endl;
        cout<<"ithSCCtoithTopo:";
        for(int i=0; i<toponodesize; i++)
        {
            cout<<ithscc_ithtopo[i]<<" ";
        }
        cout<<endl;
        cout<<"seqnumofitscchaveNE:";
        for(int i=0; i<toponodesize; i++)
        {
            cout<<ithscchaveNE_seq[i]<<" ";
        }
        cout<<endl;
        cout<<"seqscctoithtop:"<<seq<<"   ";
        for(int i=0; i<seq; i++)
        {
            cout<<seqithscchaveNE_ithtopo[i]<<" ";
        }
        cout<<endl;
    }
}
//judge the ithscc whether accessible to the destination without other scc having ENode
void fromdestdfswithoutscchaveEN(int rdv)
{
    scc_istoanodeFromDestwithouthaveEN[rdv]=true;
    for(int i=0; i<rTopoEdge[rdv].size(); i++)
    {
        if(0<ENodesnumofithSCC[rTopoEdge[rdv][i]])
        {
            scc_istoanodeFromDestwithouthaveEN[rTopoEdge[rdv][i]]=true;
        }
        if((0==ENodesnumofithSCC[rTopoEdge[rdv][i]])&&(rTopoEdge[rdv][i]!=ithnode_ithSCC[source]))
            if(false==scc_istoanodeFromDestwithouthaveEN[rTopoEdge[rdv][i]])
                fromdestdfswithoutscchaveEN(rTopoEdge[rdv][i]);
    }
}
//judge the ithnode whether accessible to the destination without other ENode
void fromdestdfswithoutEN(int rdv)
{
    istoanodeFromDestwithouthaveEN[rdv]=true;
    for(int i=0; i<SCC_rG[rdv].size(); i++)
    {
        if(true==isENode[SCC_rG[rdv][i]])
        {
            istoanodeFromDestwithouthaveEN[SCC_rG[rdv][i]]=true;
        }
        if((false==isENode[SCC_rG[rdv][i]])&&(false==istoanodeFromDestwithouthaveEN[SCC_rG[rdv][i]]))
            if(SCC_rG[rdv][i]!=source)
                fromdestdfswithoutEN(SCC_rG[rdv][i]);
    }
}
//judge the ithnode whether accessible to the destination
void fromdestdfs(int rdv)
{
    istoanodeFromDest[rdv]=true;
    for(int i=0; i<SCC_rG[rdv].size(); i++)
    {
        if((false==istoanodeFromDest[SCC_rG[rdv][i]]))
            if(SCC_rG[rdv][i]!=source)
                fromdestdfs(SCC_rG[rdv][i]);
    }
}
void rerversdfsfromdest()
{
    if(true==initsetparavalue)
    {
        memset(scc_istoanodeFromDestwithouthaveEN,0,sizeof(scc_istoanodeFromDestwithouthaveEN));
        memset(istoanodeFromDestwithouthaveEN,0,sizeof(istoanodeFromDestwithouthaveEN));
        memset(istoanodeFromDest,0,sizeof(istoanodeFromDest));
    }
    fromdestdfswithoutscchaveEN(ithnode_ithSCC[destination]);
    // if(ithnode_ithSCC[destination]!=ithnode_ithSCC[source]){
    //     istoanodeFromDestwithoutscchaveEN[ithnode_ithSCC[source]]=false;
    //  }
    fromdestdfswithoutEN(destination);
    fromdestdfs(destination);
}

//bfs parameters
vector<int>path;
bool bfsissucceedtofindpath;

struct status
{
    int id;
    int hop;
    int costsum;

    int numbofhavedENode;
    int OutNodeisNE;
    int OutNodeisE;

    int ithSCC;
    int Enodeofithscc;
    int vsnode[MAX_V];
    int seqlocithscchavedNE;
    friend bool operator<(const struct status a,const struct status b)
    {
        if(a.costsum<b.costsum)
        {
            return true;
        }
        else
        {
            if(a.costsum>b.costsum)
            {
                return false;
            }
            else
            {
                if(a.hop<b.hop)
                {
                    return true;
                }
                else
                {
                    if(a.hop>b.hop)
                        return false;
                }
            }
        }
        return false;
    }
} fq,dp,sp;
priority_queue<status>que;
int minimumcost;
int INF=11111111;
int d[MAX_V];
bool used[MAX_V];
int firstvs[MAX_V];
vector<status>preque[MAX_EV];
int PreNode[MAX_EV];
int preNode_seq[MAX_V];
int preNodesize;

typedef pair<int,int> P;
void prepushEnode()
{
    preNodesize=0;
    PreNode[preNodesize]=source;
    preNode_seq[source]=preNodesize;
    preNodesize++;
    for(int i=0; i<ENodeSize; i++)
    {
        PreNode[preNodesize]=ENode[i];
        preNode_seq[ENode[i]]=preNodesize;
        preNodesize++;
    }
    for(int i=0; i<preNodesize; i++)
    {
        fill(d,d+Node4InitInputGragh,INF);

        fill(firstvs,firstvs+Node4InitInputGragh,-1);
        priority_queue<P,vector<P>,greater<P> >midque;
        d[PreNode[i]]=0;
        midque.push(P(0,PreNode[i]));
        while(!midque.empty())
        {
            P p=midque.top();
            midque.pop();
            int v=p.second;
            if((v==destination)||((v!=PreNode[i])&&(true==isENode[v]))||(d[v]<p.first))
                continue;
            //used[v]=true;
            int u;
            for(int j=0; j<SCC_G[v].size(); j++)
            {
                u=SCC_G[v][j];
                //if(InitGraphMatix[v][u]!=0)
                if(d[u]>(d[v]+InitGraphMatix[v][u]))
                {
                    d[u]=d[v]+InitGraphMatix[v][u];
                    //d[u]=d[v]+1;
                    firstvs[u]=v;
                    midque.push(P(d[u],u));
                }
            }
        }
        int pre,next,costsum,hop;
        //int b1,b2;
        int icci=ithnode_ithSCC[PreNode[i]];
        int iccj;
        //for(int j=0; j<Node4InitInputGragh; j++)
        int j;


        for(int l=0; l<preNodesize; l++)
        {
            if(l==0)
            {
                j=destination;
            }
            else
                j=PreNode[l];
            if(i==0&&l==0)
            {
                continue;
            }
            if((firstvs[j]!=-1)&&((true==isENode[j])||j==destination))
            {
                iccj=ithnode_ithSCC[j];

                if((icci==iccj)||((ithscchaveNE_seq[icci]+1)==ithscchaveNE_seq[iccj]))
                {
                    //cout<<"i "<<icci<<"  j "<<iccj<<endl;
                    pre=j;
                    hop=0;
                    costsum=0;
                    memset(fq.vsnode,-1,sizeof(fq.vsnode));
                    //memset(fq.bit,0,sizeof(fq.bit));
                    next=firstvs[pre];
                    fq.vsnode[pre]=next;
//                b1=pre/8;
                    //  b2=pre%8;
                    //fq.bit[b1]|=(1<<b2);
                    while(-1!=next)
                    {
                        costsum+=InitGraphMatix[next][pre];
                        pre=next;

                        hop++;
                        next=firstvs[pre];
                        fq.vsnode[pre]=next;
                        /*
                         if(next!=-1){
                             b1=pre/8;
                             b2=pre%8;
                             fq.bit[b1]|=(1<<b2);
                         }
                         */
                    }
                    fq.id=j;
                    fq.hop=hop;
                    fq.costsum=costsum;
                    preque[i].push_back(fq);
                }

            }
        }

    }
}

void showbfsunit(struct status *show)
{
    cout<<"----------------------------------"<<endl;
    cout<<"Enloc: "<<node_index[show->id]<<" weight: "<<show->costsum;
    cout<<"  hop: "<<show->hop<<endl;
    int loc=show->id;
    int pre=show->vsnode[show->id];
    while(-1!=pre)
    {
        cout<<"("<<pre<<","<<loc<<") "<<endl;
        cout<<"("<<node_index[pre]<<","<<node_index[loc]<<") "<<Index4InitGraphEdgeMatix[pre][loc]<<"  "<<InitGraphMatix[pre][loc]<<endl;
        loc=pre;
        pre=show->vsnode[loc];

    }


}


bool controltime;//=true;
int ansvs[MAX_V];
void dfsfor14or15(int start,int enodenum,int costsum)//,unsigned int bit[19])
{
    if(false==controltime)
    {
        return ;
    }
    int back_vs[MAX_V];
    //unsigned int bt[19];
    memcpy(back_vs,ansvs,Node4InitInputGragh*4);
    //memcpy(bt,bit,19*4);
    /*
    if(enodenum==ENodeSize)
    {

            int cost=getlastpath(start,vs,costsum);
            if(cost>0)
            {
                if(minimumcost>cost)
                {
                    minimumcost=cost;
                    path.clear();
                    path.push_back(destination);
                    for(int i=vs[destination]; i!=-1; i=vs[i])
                    {
                        path.push_back(i);
                    }
                    bfsissucceedtofindpath=true;
                }
                controltime=false;
            }
        return ;
    }
    */
    /*
        ftime(&endtime);
        if(17<=(endtime.time-starttime.time))
        {
            controltime=false;
            return ;
        }
    */
    //memcpy(vs,back_vs,Node4InitInputGragh*4);
    int seq=preNode_seq[start];
    bool isover;
    //int b1,b2;
    for(int i=0; i<preque[seq].size(); i++)
    {
        if(true==controltime)
        {
            isover=true;
            dp=preque[seq][i];
            if(dp.id!=destination)
            {
                for(int j=0; j<Node4InitInputGragh; j++)
                {
                    if((-1!=ansvs[j])&&(-1!=(dp.vsnode[j])))
                    {
                        isover=false;
                    }
                }

                /*

                for(int j=0; j<19; j++)
                {
                    if(0!=(bit[j]&dp.bit[j]))
                    {
                        //cout<<bit[j]<<"     "<<dp.bit[j]<<"jj"<<j<<endl;
                        isover=false;
                    }
                }
                */
                if(isover)
                {
                    if((enodenum+1)==ENodeSize)
                    {
                        if(istoanodeFromDestwithouthaveEN[dp.id]&&scc_istoanodeFromDestwithouthaveEN[ithnode_ithSCC[dp.id]])
                        {
                            int start=dp.id;
                            int seqstart=preNode_seq[start];
                            int cost=costsum+dp.costsum;
                            for(int j=0; j<preque[seqstart].size(); j++)
                            {
                                sp=preque[seqstart][j];
                                if(sp.id==destination)
                                {
                                    break;
                                }

                            }

                            if(sp.vsnode[destination]!=-1)
                            {
                                for(int j=0; j<Node4InitInputGragh; j++)
                                {
                                    if(-1!=(dp.vsnode[j]))
                                    {
                                        ansvs[j]=dp.vsnode[j];
                                        //b1=j/8;
                                        // b2=j%8;
                                        //bit[b1]|=(1<<b2);
                                    }

                                }
                                controltime=false;
                                if(-1!=(sp.vsnode[destination]))
                                {
                                    int pre=destination;
                                    int next;//=firstvs[pre];
                                    do
                                    {
                                        next=sp.vsnode[pre];
                                        ansvs[pre]=next;
                                        costsum+=InitGraphMatix[next][pre];
                                        pre=next;
                                        //slpsp->hop++;
                                    }
                                    while(next!=start);
                                    if(minimumcost>cost)
                                    {
                                        minimumcost=cost;
                                        path.clear();
                                        path.push_back(destination);
                                        for(int i=ansvs[destination]; i!=-1; i=ansvs[i])
                                        {
                                            path.push_back(i);
                                        }
                                        bfsissucceedtofindpath=true;
                                    }
                                }
                            }
                        }
                    }
                    else
                    {

                        for(int j=0; j<Node4InitInputGragh; j++)
                        {
                            if(-1!=(dp.vsnode[j]))
                            {
                                ansvs[j]=dp.vsnode[j];
                                //b1=j/8;
                                // b2=j%8;
                                //bit[b1]|=(1<<b2);
                            }
                        }
                        dfsfor14or15(dp.id,enodenum+1,costsum+dp.costsum);//,bit);
                        memcpy(ansvs,back_vs,Node4InitInputGragh*4);
                        //memcpy(bit,bt,19*4);

                    }

                }

            }
        }
    }
    return ;
}
void bfsfor14or15()
{
    minimumcost=INF;
    prepushEnode();

    for(int i=0; i<preNodesize; i++)
        sort(preque[i].begin(),preque[i].end());
    if(debug_showbfsfor14or15)
        for(int i=0; i<preNodesize; i++)
        {

            cout<<"+++++++++++"<<node_index[PreNode[i]]<<"+++++++++++++"<<endl;
            for(int j=0; j<preque[i].size(); j++)
            {
                sp=preque[i][j];
                showbfsunit(&sp);
            }

        }

    // unsigned int bit[19];
    fill(ansvs,ansvs+Node4InitInputGragh,-1);
    //fill(bit,bit+19,0);
    controltime=true;
    //unsigned int bit1,bit2;
    //bit1=source/8;
    //bit2=source%8;
    //bit[bit1]|=(1<<bit2);
    dfsfor14or15(source,0,0);//,bit);
    //cout<<"minimumcost:    "<<minimumcost<<endl;

}





/************************************************************************/




void search_route(char *topo[5000], int edge_num, char *demand)
{
    int edgeNum = edge_num;
    int nodeNum  = Graph::readNodeNum(topo, edge_num);
    int includingSetSize = Graph::readConditionSize(demand);

    cout << "edgeNum = " << edgeNum << " " << "nodeNum = " << nodeNum << " " << "includingSetSize = " << includingSetSize << endl;
    //if(edgeNum >= 3000 && nodeNum >= BigNode && includingSetSize >= 30)
    if( nodeNum > BigNode)
    {
        unsigned short result[MAX_V];
        int demandStrLength = strlen(demand);
        int mInput;
        int ithtopolebgth;
        int mFlag,inEdge,outEdge,mWeight;
        // 14, 15 题代码
        /********************************************/
        ftime(&starttime);
        ENodeSize=0;
        if(debug_isruningsearch_route)
        {
            debug_showinputstr=false;
            debug_write=false;
            debug_showbfsqueue=false;
            debug_showinitmatrix=false;
            debug_showresultpath=false;
            debug_showtopoinfo=false;
            debug_showbfsfor14or15=false;
            initcleardata();

        }
        beabletofindpath=true;
        if(debug_showinputstr)
            cout<<"demand:"<<demand<<endl;
        int demandi=0;
        for(mInput=0; demandi<demandStrLength; demandi++)
        {
            if(!((demand[demandi]<='9')&&(demand[demandi]>='0')))
            {
                demandi++;
                break;
            }
            mInput=mInput*10+(demand[demandi]-'0');
        }
        source = mInput;

        for(mInput=0; demandi<demandStrLength; demandi++)
        {
            if(!((demand[demandi]<='9')&&(demand[demandi]>='0')))
            {
                demandi++;
                break;
            }
            mInput=mInput*10+(demand[demandi]-'0');
        }
        destination = mInput;

        index_node[source]=Node4InitInputGragh;
        node_index[Node4InitInputGragh]=source;
        source=Node4InitInputGragh;

        ExistNodeafterscaledown[Node4InitInputGragh]=1;
        Node4InitInputGragh++;

        index_node[destination]=Node4InitInputGragh;
        node_index[Node4InitInputGragh]=destination;
        destination=Node4InitInputGragh;

        ExistNodeafterscaledown[Node4InitInputGragh]=1;
        Node4InitInputGragh++;

        for(mInput=0; demandi<=demandStrLength; demandi++)
        {
            if(demand[demandi]!='\n')
            {
                if((demand[demandi]==',')||(demand[demandi]=='|')||(demandi==(demandStrLength)))
                {
                    index_node[mInput]=Node4InitInputGragh;
                    node_index[Node4InitInputGragh]=mInput;
                    ENode[ENodeSize]=Node4InitInputGragh;
                    isENode[Node4InitInputGragh]=true;
                    ExistNodeafterscaledown[Node4InitInputGragh]=1;
                    Node4InitInputGragh++;
                    ENodeSize++;
                    mInput=0;
                    if(demandi==demandStrLength)
                    {
                        break;
                    }
                }
                else
                {
                    mInput=mInput*10+(demand[demandi]-'0');
                }
            }
            if((demand[demandi]!='\n')&&(demandi==demandStrLength))
                break;
        }

        for(int i=0; i<edge_num; i++)
        {
            ithtopolebgth=strlen(topo[i]);
            int j=0;
            mInput=0;
            for(; j<ithtopolebgth; j++)
            {
                if(!((topo[i][j]<='9')&&(topo[i][j]>='0')))
                {
                    j++;
                    break;
                }
                mInput=mInput*10+(topo[i][j]-'0');
            }
            mFlag=mInput;
            mInput=0;
            for(; j<ithtopolebgth; j++)
            {
                if(!((topo[i][j]<='9')&&(topo[i][j]>='0')))
                {
                    j++;
                    break;
                }
                mInput=mInput*10+(topo[i][j]-'0');
            }
            if(index_node[mInput]==-1)
            {
                index_node[mInput]=Node4InitInputGragh;
                node_index[Node4InitInputGragh]=mInput;
                ExistNodeafterscaledown[Node4InitInputGragh]=1;
                Node4InitInputGragh++;
            }
            inEdge=index_node[mInput];

            mInput=0;
            for(; j<ithtopolebgth; j++)
            {
                if(!((topo[i][j]<='9')&&(topo[i][j]>='0')))
                {
                    j++;
                    break;
                }
                mInput=mInput*10+(topo[i][j]-'0');
            }

            if(index_node[mInput]==-1)
            {
                index_node[mInput]=Node4InitInputGragh;
                node_index[Node4InitInputGragh]=mInput;
                ExistNodeafterscaledown[Node4InitInputGragh]=1;
                Node4InitInputGragh++;

            }
            outEdge=index_node[mInput];
            mInput=0;
            for(; j<ithtopolebgth; j++)
            {
                if(!((topo[i][j]<='9')&&(topo[i][j]>='0')))
                {
                    j++;
                    break;
                }
                mInput=mInput*10+(topo[i][j]-'0');
            }
            mWeight=mInput;
            if((inEdge!=destination)&&(outEdge!=source)&&(inEdge!=outEdge))
            {

                if(InitGraphMatix[inEdge][outEdge]==0)
                {
                    InitGraphMatix[inEdge][outEdge]=mWeight;
                    Index4InitGraphEdgeMatix[inEdge][outEdge]=mFlag;
                    Edg4InitInputGragh++;

                }
                else
                {
                    if(InitGraphMatix[inEdge][outEdge]>mWeight)
                    {
                        InitGraphMatix[inEdge][outEdge]=mWeight;
                        Index4InitGraphEdgeMatix[inEdge][outEdge]=mFlag;
                    }
                }

            }
        }

        scaledowngragh();
        if(beabletofindpath==true)
            beabletofindpath=judgegraph();

        convertMatrixtoList();

        SCCSize=scc();

        if(debug_showinputstr)
            DebugPrint();
        if((isDFStodestination==false)||(SCCDFSAllENode!=ENodeSize))
            beabletofindpath=false;

        //after SCC,the graph's node accimulate to some parts,some parts compose
        //the DAG,get the toposort.
        getTopoSort4SCC();

        rerversdfsfromdest();
        bfsfor14or15();
        beabletofindpath=bfsissucceedtofindpath;
        int pathcostsum=0;
        if(debug_showresultpath)
        {
            for(int i=path.size()-1; i>=1; i--)
            {
                if(i==(path.size()-1))
                {
                    cout<<"edge index:"<<Index4InitGraphEdgeMatix[path[i]][path[i-1]];
                    pathcostsum+=InitGraphMatix[path[i]][path[i-1]];
                }
                else
                {
                    cout<<"|"<<Index4InitGraphEdgeMatix[path[i]][path[i-1]];
                    pathcostsum+=InitGraphMatix[path[i]][path[i-1]];
                }
            }
            cout<<endl<<"pathcostsum:"<<pathcostsum<<endl;
        }
        if(debug_showresultpath)
        {
            cout<<endl;
            for(int i=path.size()-1; i>=0; i--)
            {

                if(i==(path.size()-1))
                {
                    cout<<"node index:"<<node_index[path[i]];
                }
                else
                {
                    cout<<"|"<<node_index[path[i]];
                }
            }


        }
        if(true==beabletofindpath)
        {
            int j=0;
            for(int i=path.size()-1; i>=1; i--,j++)
            {
                if((1==i)&&debug_showresultpath)
                    cout<<"succeed"<<endl;
                result[j]=Index4InitGraphEdgeMatix[path[i]][path[i-1]];
            }
            for (int i = 0; i < j; i++)
                record_result(result[i]);
        }
        /**********************************************/

    }
    else
    {
        Graph linkGraph;
        InGraph inLinkGraph;
        Graph::readGraph(topo, edge_num, linkGraph, inLinkGraph);
        Graph::readCondition(demand, linkGraph);
        Graph::preprocess(linkGraph, inLinkGraph);
        vector<int> result;
        Graph::kShortestPath (result, linkGraph);

        int j = result.size();
        while(j > 0)
        {
            record_result(result[--j]);
        }
    }
}

