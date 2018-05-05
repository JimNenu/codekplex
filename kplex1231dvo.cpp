/************************************************

** This is an exact solver for maximum k-plexes in massive graphs.    
** Date:	2017.12.31  
** Author:	Jiejiang chen, Jian Gao   
                         
** To compile the solver use command: g++ -std=c++11 -O3 kplex1231dvo.cpp -o kplex
				 
************************************************/


#include <iostream>
#include <cstring>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <queue>
#include <algorithm>
#include<unordered_set>
#include<sys/time.h>
#include<sys/types.h>
#include<sys/resource.h>
#include<limits.h>
#include "Array.h"

//  #define DEBUG

using namespace std;

#define mem_2G 536870912
#define insert_v(end, value) *(end++) = value;
#define delete_i(index, end) *index = *(--end);


int v_n;//number of vertex
int e_n;//number of edge

struct Edge
{
    int v1,v2;
};

int run_time = 1000;//time_limit


struct Vertex_sort//
{
    int index;
    int degree;
};

struct Vertex//
{
    int state;
    int degree;
    int *neighbor;
};

int is_reduce = 1;
//----------
int *neighbor_len;//
int **neighbor;//
//int *state;//the sate of vertex
//---------------

int *temp_array;//
int *temp_mark;//
int *temp_index;//
int *neighbor_in_solution;//
//int *kv;
//int *v_weight;
Edge *edge;//
Vertex_sort *vertex_sort;//
Vertex *vertex;//

Array *crusol;//
Array *remaining_vertex;//
int para_k = 1;//

Array *U;//
Array *S;//  

int *best_solution;//
int best_solution_size = 0;//

int max_degree = -1;//

double read_time, init_time, search_time;//
int UB = INT_MAX;//
int LB = 0;//
int *t_U;//

int **block;//
int cur_block = -1;//
int block_size = 0;//

unordered_set<unsigned long long> edge_hash;//

void freeAll();//

bool cmp (Vertex_sort v1, Vertex_sort v2)//
{
    return v1.degree > v2.degree;
}

double get_utime()//
{
    struct rusage utime;
    getrusage(RUSAGE_SELF, &utime);
    return (double) (utime.ru_utime.tv_sec + (double)utime.ru_utime.tv_usec / 1000000);
}

unsigned long long encode_pairID(unsigned long long  v1, unsigned long long  v2)//
{
    unsigned long long n1,n2;
    unsigned long long pairID;
    if(v1 < v2)
    {
        n1 = v1; n2 = v2;
    }
    else
    {
        n1 = v2; n2 = v1;
    }
    pairID = ((n1 + n2 + 1) * (n1 + n2) >> 1) + n2;
    return pairID;
}

int edge_is(int m, int n)//
{
    unsigned long long id = encode_pairID(m, n);
    if(edge_hash.count(id))
        return 1;
    else
        return 0;
}

int* expand_memory()//
{
    if(cur_block == block_size - 1)
    {
        block[block_size] = (int *)malloc(mem_2G * sizeof(int));
        if(block[block_size] == NULL)
        {
            freeAll();
            exit(0);
        }
        block_size++;
    }
    cur_block++;
    return  block[cur_block];
}

void readGraph(char* File_Name)//
{

    ifstream FIC;
    FIC.open(File_Name);
    if(FIC.fail())
    {
        printf("### Error open, File_Name %s\n", File_Name);
    }
    char strReading[1024];
    char p[10], temp[10];
    FIC >> p >> temp >> v_n >> e_n;

    neighbor = (int **)malloc(v_n * sizeof(int**));
    block = (int **)malloc(100 * sizeof(int**));
    neighbor_len = (int *)malloc(v_n * sizeof(int));
    vertex_sort= (Vertex_sort *)malloc(v_n *sizeof(Vertex_sort));
    edge = (Edge *)malloc(e_n * sizeof(Edge));
    temp_array = (int *)malloc(v_n * sizeof(int));
    temp_mark = (int *)malloc(v_n * sizeof(int));
    temp_index = (int *)malloc(v_n * sizeof(int));
    vertex = (Vertex *)malloc(v_n * sizeof(Vertex));
    neighbor_in_solution = (int *)malloc(v_n * sizeof(int));
    t_U = expand_memory();
    best_solution = (int *)malloc(v_n *sizeof(int));

    remaining_vertex = new Array(v_n);///////////////////////////////////////////
    crusol = new Array(v_n);
    U = new Array(v_n);

    for(int i = 0; i < v_n; i++)
    {
         vertex[i].degree = 0;
         vertex[i].state = 0;
    }

    int a,b;
    int count = 0;
    unsigned long long id;

    while(FIC >> temp >> a >> b)
    {
        if(strcmp(temp, "e") ==  0)
        {
            a--;
            b--;
            vertex[a].degree++;
            vertex[b].degree++;
            edge[count].v1 = a;
            edge[count].v2 = b;
            id = encode_pairID(a, b);
            edge_hash.insert(id);
            count++;
        }
    }

    for(int i = 0; i < v_n; i++)
    {
        if(vertex[i].degree > max_degree)
            max_degree = vertex[i].degree;
        vertex[i].neighbor = (int *)malloc(vertex[i].degree * sizeof(int));
        neighbor[i] = (int *)malloc(vertex[i].degree * sizeof(int));
        remaining_vertex->insert_element(i);
    }

    for(int i = 0; i < v_n; i++)
        vertex[i].degree = 0;
    for(int e  = 0; e < e_n; e++)
    {
        a = edge[e].v1;
        b = edge[e].v2;
        vertex[a].neighbor[vertex[a].degree] = b;
        vertex[b].neighbor[vertex[b].degree] = a;
        vertex[a].degree++;
        vertex[b].degree++;
    }
    free(edge);
    FIC.close();
    printf("### the number of vertex:\t%d\n### the number of edge:\t\t%d\n### maximum degree:\t\t%d\n### time limit: \t\t%ds\n",
           v_n, e_n, max_degree, run_time);

}

void freeAll()//
{
    free(neighbor_len);
    free(temp_array);
    free(temp_mark);
    free(temp_index);
    free(neighbor_in_solution);
    free(best_solution);
    for(int i = 0; i <100; i++)
        free(block[i]);
    free(block);
    for(int i = 0; i < v_n; i++)
    {
        free(vertex[i].neighbor);
        free(neighbor[i]);
    }
    free(vertex);
    free(neighbor);
    free(vertex_sort);
}

void update_best_solution()//
{
    best_solution_size = crusol->size();
    for(int i = crusol->begin(); i < crusol->size(); i++)
        best_solution[i] = crusol->element_at(i);
    LB = best_solution_size;
}

bool check_solution()//
{
    int i, j, adj_num, v, u;
    for(i = 0; i < best_solution_size; i++)
    {
        adj_num = 0;
        v = best_solution[i];
        for(j = 0; j < best_solution_size; j++)
        {
            u = best_solution[j];
            if(v ==u)
                continue;
            if(edge_is(v,u))
                adj_num++;
        }
        if(adj_num < best_solution_size - para_k)
            return false;
    }
    return true;
}

void printf_solution()//
{
    if(check_solution())
    {
        printf("after checking , the solution is correct,  solution size: %d, time: %f\n", best_solution_size, get_utime());
    }
    else
    {
        printf("the solution found is wrong\n");
    }

}

/*void reduce_edge_from_neighbor(int v, int mark)//update the neighbor list//
{
    int i,j;
    int u;
    if(mark == 1)
    {
        for(i = 0; i < vertex[v].degree; i++)
        {
            u = vertex[v].neighbor[i];
            for(j = 0; j < vertex[u].degree; j++)
                if(vertex[u].neighbor[j] == v)
                    break;
            vertex[u].degree--;
            vertex[u].neighbor[j] = vertex[u].neighbor[vertex[u].degree];
        }
        vertex[v].degree = 0;
    }
    else if(mark == 2)
    {
        for(i = 0; i < neighbor_len[v]; i++)
        {
            u = neighbor[v][i];
            for(j = 0; j < neighbor_len[u]; j++)
                if(neighbor[u][j] == v)
                    break;
            neighbor_len[u]--;
            neighbor[u][j] = neighbor[u][neighbor_len[u]];
        }
        //neighbor_len[v] = -1;
    }
}*/

void reduce_graph_1()//
{
    int *degree_counter, *where, *candidate;
    int node, p1, i, j, h, k, t, neighbors,_tsize;
    int u,v;
    int cur_degree = 0;
    degree_counter = temp_mark;
    where = temp_index;
    candidate = temp_array;

    for(i = 0; i <= max_degree; i++)
        degree_counter[i] = 0;
    for(node = remaining_vertex->begin(); node < remaining_vertex->size(); node++)
    {
        v = remaining_vertex->element_at(node);
        degree_counter[vertex[v].degree]++;
    }
    j = 0;
    for(i = 0; i <= max_degree; i++)
    {
        k = degree_counter[i];
        degree_counter[i] = j;
        j += k;
    }

    for(node = remaining_vertex->begin(); node < remaining_vertex->size(); node++)
    {
        v = remaining_vertex->element_at(node);
        t = degree_counter[vertex[v].degree];
        degree_counter[vertex[v].degree]++;
        candidate[t] = v;
        where[v]  = t;
    }

    for(i = max_degree; i > 0; i--)
        degree_counter[i] = degree_counter[i - 1];

    degree_counter[0] = 0;
    p1 = 0;
    cur_degree = vertex[candidate[0]].degree;

    for(i = remaining_vertex->begin(); i < remaining_vertex->size(); i++)
    {
        v = remaining_vertex->element_at(i);
        neighbor_len[v] = vertex[v].degree;
    }

    while(p1 < remaining_vertex->size())
    {
        node = candidate[p1];
        if(p1 < remaining_vertex->size() - 1  && neighbor_len[node] == neighbor_len[candidate[p1+1]])
            degree_counter[neighbor_len[node]] = p1 + 1;
        if(neighbor_len[node] + para_k - 1 >= remaining_vertex->size() - p1 - 1)
        {
            crusol->clear();
            for(i = p1; i < remaining_vertex->size(); i++)
                crusol->insert_element(candidate[i]);
            break;
        }

        for(i = 0; i < vertex[node].degree; i++)
        {
            neighbors = vertex[node].neighbor[i];
            if(where[neighbors] > p1)
            {
                t = where[neighbors];
                h = degree_counter[neighbor_len[neighbors]];
                k = candidate[h];

                candidate[h] = neighbors;
                where[neighbors] = h;
                candidate[t] = k;
                where[k]= t;

                degree_counter[neighbor_len[neighbors]]++;
                neighbor_len[neighbors]--;
                if(neighbor_len[neighbors] != neighbor_len[candidate[h-1]])
                    degree_counter[neighbor_len[neighbors]] = h;
            }
        }
        p1++;
    }

    if(crusol->size() > best_solution_size)
        update_best_solution();

    queue<int> remove_que;
    for(i = remaining_vertex->begin(); i < remaining_vertex->size(); i++)
    {
        v = remaining_vertex->element_at(i);
        neighbor_len[v] = vertex[v].degree;
        if(vertex[v].degree + para_k <= best_solution_size ||  remaining_vertex->size() <= best_solution_size)
        {
            vertex[v].state = -1;
            remove_que.push(v);
        }
    }

    while(!remove_que.empty())//
    {
            v = remove_que.front();
            remove_que.pop();
            remaining_vertex->delete_element(v);
            for(i = 0; i < vertex[v].degree; i++)
            {
                u = vertex[v].neighbor[i];
                neighbor_len[u]--;
                if(!vertex[u].state && (neighbor_len[u] + para_k <= best_solution_size ||  remaining_vertex->size() <= best_solution_size))
                {
                    vertex[u].state = -1;
                    remove_que.push(u);
                }
            }
    }

    int mark = 1;
    for(i = remaining_vertex->begin(); i < remaining_vertex->size(); i++)//
    {
        v = remaining_vertex->element_at(i);
        mark = 1;
        for(j = 0; j < vertex[v].degree; j++)
        {
            u = vertex[v].neighbor[j];
            if(vertex[u].state)
            {
                mark = 0;
                vertex[v].degree--;
                vertex[v].neighbor[j] = vertex[v].neighbor[vertex[v].degree];
                j--;
            }
        }
        if(!mark && U->is_in_array(v))//
            U->delete_element(v);
    }
}


void  reduce_graph_2(Array *Uset)//
{
    int i,v,j,k,p,u1,u2,tnode;
    int rv, rd, cu, kv, ku;
    int ssize;
    ssize = 0;//
    queue<int> remove_que;
    for(i = Uset->begin(); i < Uset->size(); i++)
    {
        u1 = Uset->element_at(i);
        neighbor_len[u1] = 0;
    }

    for(i = Uset->begin(); i < Uset->size(); i++)
    {
        v = Uset->element_at(i);//
        if(U->is_in_array(v))//
            continue;
        kv = para_k - 1;//
        rv = LB - kv - ssize;
        rd = vertex[v].degree;
        neighbor_len[v] = 0;
        for(j = 0; j < rd; j++)//
        {
            u1 = vertex[v].neighbor[j];
            if(!Uset->is_in_array(u1))
                continue;
            neighbor_len[u1] = 1;
            neighbor[u1][0] = v;
            //neighbor_len[v]++;
            neighbor[v][neighbor_len[v]++] = u1;
        }
        if(neighbor_len[v]  < rv)//
        {
            Uset->delete_element(v);
            //reduce_edge_from_neighbor(v, 1);
            i--;
            continue;
        }
        while(!remove_que.empty())//
            remove_que.pop();
        for(j = 0; j < rd; j++)//
        {
            u1 = vertex[v].neighbor[j];
            if(!Uset->is_in_array(u1))//
                continue;
            for(k = j + 1; k < rd; k++)
            {
                u2 = vertex[v].neighbor[k];
                if(!Uset->is_in_array(u2))
                    continue;
                if(edge_is(u1, u2))
                {
                    neighbor[u1][neighbor_len[u1]++] = u2;
                    neighbor[u2][neighbor_len[u2]++] = u1;
                }
            }
            ku = para_k - 1;
            cu = neighbor_len[u1] + ku < neighbor_len[v]? neighbor_len[u1] + ku : neighbor_len[v];
            temp_index[u1] = neighbor_len[u1];//
            if(cu < rv)//
            {
                temp_index[u1] = 0;//
                remove_que.push(u1);
            }
        }

        temp_index[v] = neighbor_len[v];//
        while(!remove_que.empty())//
        {
            u1 = remove_que.front();
            remove_que.pop();
            for(j = 0; j < neighbor_len[u1]; j++)
            {
                u2 = neighbor[u1][j];
                if(!temp_index[u2])//
                    continue;
                temp_index[u2]--;//
                if(temp_index[u2] + ku < rv && u2 != v)//
                {
                    temp_index[u2] = 0;//
                    remove_que.push(u2);
                }
            }
            if(temp_index[v] < rv)//
            {
                Uset->delete_element(v);
                i--;
                break;
            }
        }
    }

    U->clear();//
    int mark = 1;
    for(i = Uset->begin(); i < Uset->size(); i++)//
    {
        v = Uset->element_at(i);
        mark = 1;
        for(j = 0; j < vertex[v].degree; j++)
        {
            u1 = vertex[v].neighbor[j];
            if(!Uset->is_in_array(u1))
            {
                mark = 0;
                vertex[v].degree--;
                vertex[v].neighbor[j] = vertex[v].neighbor[vertex[v].degree];
                j--;
            }
        }
        if(mark)//
            U->insert_element(v);
    }
}

void  reduce_graph_in_BB(int* &begin, int* &end, int rev)//
{
    int i,v,j,k,p,u1,u2,tnode;
    int rv, rd, cu, kv, ku;
    int mark = 1;
    int ssize;
    int *ii;
    int t_LB = -1;
    queue<int> remove_que;
    ssize = S->size();//

    for(ii = begin; ii < end; ii++)
        temp_index[*ii] = 1;//

    for(i = 0; i < vertex[rev].degree; i++)
    {
        u1 = vertex[rev].neighbor[i];
        if(temp_index[u1])//
            temp_index[u1] = 2;
    }

    for(ii = begin; ii < end; ii++)
    {
        mark = 1;
        v  = *ii;
        if(temp_index[v] < 2)
            continue;
        kv = para_k - 1 - (ssize - neighbor_in_solution[v]);//
        /*if(kv < 0)//
        {
            cout << "@";
            delete_i(ii, end);
            ii--;
            temp_index[v] = 0;//
            /*for(i = 0; i < vertex[v].degree; i++)
            {
                u1 = vertex[v].neighbor[i];
                if(temp_index[u1])
                    vertex[u1].state--;
            }
            continue;
        }*/
        rv = LB - kv - ssize;//
        rd = vertex[v].degree;
        neighbor_len[v] = 0;
        for(j = 0; j < rd; j++)//
        {
            u1 = vertex[v].neighbor[j];
            if(!temp_index[u1])
                continue;
            neighbor[u1][0] = v;
            neighbor_len[u1] = 1;
            neighbor[v][neighbor_len[v]++] = u1;
        }

        if(neighbor_len[v]  < rv)//
        {
            delete_i(ii, end);
            ii--;
            temp_index[v] = 0;
          /*  for(i  = 0; i < neighbor_len[v]; i++)
            {
                u1 = neighbor[v][i];
                vertex[u1].state--;
            }*/
            continue;
        }

        while(!remove_que.empty())
            remove_que.pop();

        for(j = 0; j < rd; j++)//
        {
            u1 = vertex[v].neighbor[j];
            if(!temp_index[u1])
                continue;
            for(k = j + 1; k < rd; k++)
            {
                u2 = vertex[v].neighbor[k];
                if(!temp_index[u2])
                    continue;
                if(edge_is(u1, u2))
                {
                    neighbor[u1][neighbor_len[u1]++] = u2;
                    neighbor[u2][neighbor_len[u2]++] = u1;
                }
            }

            ku = para_k - 1 - (ssize - neighbor_in_solution[u1]);
            temp_array[u1] = neighbor_len[u1] + ku < neighbor_len[v]? neighbor_len[u1] + ku : neighbor_len[v];
            //temp_array[u1] = neighbor_len[u1] + ku;
            if(temp_array[u1] < rv )
            {
                temp_array[u1] = 0;//
                remove_que.push(u1);
            }
        }

        temp_array[v] = neighbor_len[v];//
        while(!remove_que.empty())//
        {
            u1 = remove_que.front();
            remove_que.pop();
            temp_array[v]--;//
            for(j = 1; j < neighbor_len[u1]; j++)//
            {
                u2 = neighbor[u1][j];
                if(!temp_array[u2])//
                    continue;
                temp_array[u2]--;//
                if(temp_array[u2] > temp_array[v])//
                    temp_array[u2] = temp_array[v];
                if(temp_array[u2] < rv)//
                {
                    temp_array[u2] =0;
                    remove_que.push(u2);
                }
            }
            if(temp_array[v] < rv)//
            {
                delete_i(ii, end);
                ii--;//
                temp_index[v] = 0;//
            /*  for(i = 0; i < neighbor_len[v]; i++)
                {
                    u1 = neighbor[v][i];
                    vertex[u1].state--;
                }*/
                //mark = 0;
                break;
            }
        }
    }

    for(ii = begin; ii < end; ii++)//
    {
        v = *ii;
        vertex[v].state = 0;
        for(i = 0; i < vertex[v].degree; i++)
        {
            if(temp_index[vertex[v].neighbor[i]])
                vertex[v].state++;
        }
    }
}

void reduce_graph()//
{
    int rn = remaining_vertex->size();
    printf("best solution\tcurrent solution\tremaining vertex\n");
    printf("%8d\t%8d\t\t%8d\n", 0, 0, v_n);
    while(true)//
    {
        reduce_graph_1();
        printf("%8d\t%8d\t\t%8d(reduce1)\n", best_solution_size, crusol->size(), remaining_vertex->size());
        reduce_graph_2(remaining_vertex);
        printf("%8d\t%8d\t\t%8d(reduce2)\n", best_solution_size, crusol->size(), remaining_vertex->size());

        if(remaining_vertex->size() < rn)
            rn = remaining_vertex->size();
        else
            break;
    }

}

int calculate_upbound(Array *ary)//
{

    int i,vd,v,j;
    for(i = ary->begin(); i < ary->size(); i++)
    {
        vd = 0;
        v = ary->element_at(i);
        vertex_sort[i].index = v;
        vertex_sort[i].degree = vertex[v].degree;
    }
    sort(vertex_sort, vertex_sort + ary->size(), cmp);
    for(i = ary->begin(); i < ary->size(); i++)
    {
        vd = vertex_sort[i].degree;
        if(i + 1 >= vd + para_k - 1)
            break;
    }
    return i + 1;
}

int calculate_upbound_in_BB(int *begin, int *end)//
{
    int i = 0,vd,v,j;
    int *ii;
    int tt = 0;
    for(ii = begin; ii < end; ii++)
    {
        v = *ii;
        //vertex_sort[i].index = v;
        //vertex_sort[i].degree = vertex[v].state + para_k  - (S->size() - neighbor_in_solution[v]);//
        //vertex_sort[i].degree = vd + para_k  - (S->size() - neighbor_in_solution[v]);
        neighbor_len[v] = vertex[v].state;
        if(neighbor_len[v] + para_k + neighbor_in_solution[v] > LB)
            tt++;
        temp_index[v] = 0;
    }
    return tt;

    //sort(vertex_sort, vertex_sort + i, cmp);//
    //for(ii = begin; ii < end; ii++)
        //temp_index[*ii] = 0;
   // for(i = 0; i < remaining_vertex->size(); i++)
       // if(temp_index[remaining_vertex->element_at(i)])
          //  cout << "$";
    /*for(j = 0; j < i; j++)
    {
        vd = vertex_sort[j].degree;
        if(j + 1 >= vd)
            break;
    }*/

   // if(j + 1 < LB && end-begin > LB)
    //cout << j+1 << "&" << end - begin << "&" << tt <<"&" << LB - S->size() << " ";
    //return j + 1;
    //return tt;

}

void reduceU1(int* &tbegin, int* &tend)
{
    int *ii, v, i, u;
    queue<int> remove_que;
    for(ii = tbegin; ii < tend; ii++)//
        temp_index[*ii] = 1;

    for(ii = tbegin; ii < tend; ii++)
    {
        v = *ii;
        vertex[v].state = 0;
        for(i = 0; i < vertex[v].degree; i++)
            if(temp_index[vertex[v].neighbor[i]] > 0)
                vertex[v].state++;//
        if(vertex[v].state + neighbor_in_solution[v] + para_k <= LB)//
        {
            remove_que.push(v);
            temp_index[v] = 2;
        }
    }

    while(!remove_que.empty())//
    {
        v = remove_que.front();
        remove_que.pop();
        temp_index[v] = 0;
        for(i = 0; i < vertex[v].degree; i++)
        {
            u = vertex[v].neighbor[i];
            if(temp_index[u] != 1)
                continue;
            vertex[u].state--;
            if(vertex[u].state + neighbor_in_solution[u] + para_k <= LB)
            {
                temp_index[u] = 0;
                remove_que.push(u);
            }
        }
    }

    for(ii = tbegin; ii < tend; ii++)//
    {
        v = *ii;
        if(!temp_index[v])
        {
            delete_i(ii, tend);
            ii--;
        }
    }
}

void reduceV(int* &begin, int* &end, int v)//
{
    int i,j,k,p,u1,u2,tnode;
    int rv, rd, cu, kv, ku;
    int ssize;
    int *ii;
    queue<int> remove_que;
    ssize = S->size();
    kv = para_k - (ssize - neighbor_in_solution[v]);//
    rv = LB - kv - ssize;//
    rd = vertex[v].degree;
    neighbor_len[v] = 0;
    for(j = 0; j < rd; j++)//
    {
        u1 = vertex[v].neighbor[j];
        if(!temp_index[u1])
            continue;
        neighbor[u1][0] = v;
        neighbor_len[u1] = 1;
        neighbor[v][neighbor_len[v]++] = u1;
    }
    for(j = 0; j < rd; j++)//
    {
        u1 = vertex[v].neighbor[j];
        if(!temp_index[u1])
            continue;
        for(k = j + 1; k < rd; k++)
        {
            u2 = vertex[v].neighbor[k];
            if(!temp_index[u2])
                continue;
            if(edge_is(u1, u2))
            {
                neighbor[u1][neighbor_len[u1]++] = u2;
                neighbor[u2][neighbor_len[u2]++] = u1;
            }
        }
        ku = para_k - 1 - (ssize - neighbor_in_solution[u1]);
        temp_array[u1] = neighbor_len[u1] + ku < neighbor_len[v]? neighbor_len[u1] + ku : neighbor_len[v];
        if(temp_array[u1] < rv)//
        {
            temp_array[u1] = 0;
            remove_que.push(u1);
            //cout << "%";
        }
    }
    temp_array[v] = neighbor_len[v];

    while(!remove_que.empty())//
    {
        u1 = remove_que.front();
        remove_que.pop();
        temp_array[v]--;//
        temp_index[u1] = 0;
        for(j = 0; j < vertex[u1].degree; j++)
        {
            u2 = vertex[u1].neighbor[j];//
            if(temp_index[u2])
                vertex[u2].state--;
        }
        for(j = 1; j < neighbor_len[u1]; j++)//
        {
            u2 = neighbor[u1][j];
            //vertex[u2].state--;
            if(!temp_array[u2])
                continue;
            temp_array[u2]--;
            if(temp_array[u2] > temp_array[v])
                temp_array[u2] = temp_array[v];
            if(temp_array[u2] < rv)
            {
                temp_array[u2] = 0;
                remove_que.push(u2);
            }
        }
    }

    for(ii = begin; ii < end; ii++)//
    {
        u1 = *ii;
        if(!temp_index[u1])
        {
            delete_i(ii, end);
            ii--;
        }
    }
}

int is_S;
int nn= 0;
void update_candidate(int* &tbegin, int* &tend, int mark, int rev)//
{
    int tt1,tt2 = 0,tt3 = 0,tt4;
    int v, i, j;
    int *ii;
    //tt1 = tend-tbegin;
    if(mark)//
    {
        //cout << "\npop:";
        reduceU1(tbegin, tend);
        //tt2 = tend-tbegin;
        reduce_graph_in_BB(tbegin, tend, rev);
    }
    else//
    {
        //cout << "\npush:";
        for(ii = tbegin; ii < tend; ii++)
            temp_mark[*ii] = 0;
        is_S = 0;
        for(i = S->begin(); i < S->size(); i++)//
        {
            v = S->element_at(i);
            if(neighbor_in_solution[v] + para_k == S->size())
            {
                for(j = 0; j < vertex[v].degree; j++)
                        temp_mark[vertex[v].neighbor[j]]++;
                is_S++;
            }
        }
        for(ii = tbegin; ii < tend; ii++)
        {
            v = *ii;
            if((is_S && temp_mark[v] != is_S) || neighbor_in_solution[v] + para_k <= S->size())
            {
                delete_i(ii, tend);
                //temp_index[v] = 0;
                ii--;
            }
        }

        //tt2 = tend - tbegin;
        reduceU1(tbegin, tend);
       //tt3 = tend - tbegin;
        reduceV(tbegin, tend, rev);
    }
    /*if(mark)
    {
        reduce_graph_in_BB(tbegin, tend, rev);
    }
    else
    {
        reduceV(tbegin, tend, rev);
    }*/
    //tt4 = tend - tbegin;
    //if(mark == 0 && tt2 == tt1 && tt3 < tt2)
        //cout << tt1 << "@" << tt2 << "@" << tt3 << "@" << tt4 << " ";
}

int ct =  0;//

int BB(int* &tbegin, int* &tend, int *tv, int nt)//
{

    if(get_utime() > run_time)//
    {
        freeAll();
        printf("out of time\n");
        exit(0);
    }
    int u, i, upper, lower, maxd, v, vd, d;
    int *maxv;
    int *ii;
    int *begin, *end;
    int *cur_ind = block[cur_block];//
    int t_cur_block = cur_block;
    if(tbegin < tend)
    {
        ct++;
        nt++;
        v = *tv;
        delete_i(tv, tend);//
        S->insert_element(v);//
        for(i = 0; i < vertex[v].degree; i++)//
            neighbor_in_solution[vertex[v].neighbor[i]]++;
        begin = tend + 1;//
        if(tend - tbegin > cur_ind + mem_2G - tend - 2)//
        {
            begin = expand_memory();
        }
        end = begin;
        for(ii = tbegin; ii < tend; ii++)//
        {
            insert_v(end, *ii);
        }

       // cout << end - begin << "&";
        //reduceV(begin, end, v);
        //cout << end - begin << " ";
        update_candidate(begin, end, 0, v);//
        upper = calculate_upbound_in_BB(begin, end) + S->size();//
        if(upper > LB)
        {
            maxd = -1;
			int minkv=para_k;
            for(ii = begin; ii < end; ii++)//
            {
                u = *ii;
                //d = neighbor_len[u] + para_k - 1 - (S->size() - neighbor_in_solution[u]);//
                d = neighbor_len[u];//
                int kkk=para_k - 1 - (S->size() - neighbor_in_solution[u]);
				if(kkk==0 )
				{
					if(minkv==0 && d>maxd || minkv>0)
					{
					   maxv = ii;
					   minkv= kkk;
					   maxd = d;
					}

				}

				else if(minkv>0 && d > maxd)
                {

                    maxd = d;
                    maxv = ii;
					minkv=kkk;
                }

				/*if(d > maxd)
				{
					maxd = d;
                    maxv = ii;

				}*/

            }

            lower = BB(begin, end, maxv, nt);
            if(lower > LB)//
            {
                LB = lower;
                update_best_solution();
                printf_solution();
            }
        }
        S->delete_element(v);//

        for(i = 0; i < vertex[v].degree; i++)//
            neighbor_in_solution[vertex[v].neighbor[i]]--;
        cur_block = t_cur_block;//
        begin = tend + 1;
        if(tend - tbegin > cur_ind + mem_2G - tend - 2)
        {
            begin = expand_memory();
        }
        end = begin;
        for(ii = tbegin; ii < tend; ii++)//
            insert_v(end, *ii);
        update_candidate(begin, end, 1, v);//
        upper = calculate_upbound_in_BB(begin, end) + S->size();
        if(upper > LB)
        {
            maxd = -1;
			int minkv=para_k;
            for(ii = begin; ii < end; ii++)//
            {
                u = *ii;
                //d = neighbor_len[u] + para_k - 1 - (S->size() - neighbor_in_solution[u]);//
                d = neighbor_len[u];//
                int kkk=para_k - 1 - (S->size() - neighbor_in_solution[u]);
				if(kkk==0 )
				{
					if(minkv==0 && d>maxd || minkv>0)
					{
						maxv = ii;
						minkv= kkk;
						maxd = d;
					}

				}

				else if(minkv>0 && d > maxd)
                {

                    maxd = d;
                    maxv = ii;
					minkv=kkk;
                }

				/*if(d > maxd)
				{
					maxd = d;
                    maxv = ii;

				}*/
            }
            lower = BB(begin, end, maxv, nt);
            if(lower > LB)
            {
                LB = lower;
                update_best_solution();
                printf_solution();
            }
        }
        return LB;
    }
    else
    {
        //cout << nt << " ";
        return S->size();
    }

}

void search()
{
    int i, v,  maxd = -1;
    int *begin, *end;
    int *maxv;
    end = t_U;
    begin = t_U;
    int en = 0;
    for(i = remaining_vertex->begin(); i < remaining_vertex->size(); i++)
    {
        v = remaining_vertex->element_at(i);
        //U->insert_element(v);
        neighbor_in_solution[v] = 0;
        insert_v(end, v);//
        temp_index[v] = 0;
        en += vertex[v].degree;
        if(vertex[v].degree > maxd)//
        {
            maxd = vertex[v].degree;
            maxv = end;
        }
    }
    crusol->clear();//
    S = crusol;
   // cout << "after reduce, the edge num = " << en / 2 << "\n";
    LB = BB(begin, end, maxv, 0);
}

int main(int argc, char *argv[])
{
    char* File_Name;
    File_Name = argv[1];//
    para_k = atoi(argv[2]);//
    run_time = atoi(argv[3]);//
    printf("%s: \n", File_Name);
    readGraph(File_Name);//
    read_time = get_utime();
    printf("----------------------------readGraph finish----------------------------\n");
    reduce_graph();//
    init_time = get_utime() - read_time;
    LB = best_solution_size;
    //UB = calculate_upbound(remaining_vertex);
    printf_solution();
    search();//
    printf("read time: %f; init time: %f; search time: %f\n", read_time, init_time, get_utime() - read_time - init_time);
    printf("total time: %f\n", get_utime());
    //cout << ct << endl;
   // cout << LB << endl;
    //cout << nn << endl;
    printf("------------------------------------------------------------------------\n");
    freeAll();
    return 0;
}
