#ifndef _SBUNDLE_TOOL_H
#define _SBUNDLE_TOOL_H

#include <stdio.h>
#include <utility>
#include <fstream>
#include <bits/stdc++.h>

extern unsigned int *head;
extern unsigned int *que, *nV, n;
const unsigned int INF = 0x3f3f3f3f;
using ui = unsigned int;
using ept = unsigned int;

class MaximumFlow
{
private:
    struct Node
    {
        ui from, to, next;
        ui cap;
    };
    Node *edge;
    ui *cap_backup;
    ui tol;
    ui *head;
    ui *dep;
    ui *gap, *que;
    ui *cur;
    ui *S;
    ui n;
    void BFS(ui start, ui end);
    ui *nV, *oID;
    bool flag;

public:
    MaximumFlow();
    ~MaximumFlow();
    void reserve(ui n, ui m);
    void init(ui _n);
    void addedge(ui u, ui v, ui w);
    ui SAP(ui start, ui end, ui nn);
    void BFSMinCut(ui start, std::vector<bool> &visited);
    void prepare(ui n, ui m);
    bool verify_SBundle_by_MaxFlowAlg0(ui *ids, ui R_end, char *matrix, long long row, int s, bool is_global);
    bool verify_SBundle_by_MaxFlowAlg(std::vector<ui> ids, ui R_end, char *matrix, long long row, int s, bool is_global);
    bool verify_Sbundle_by_maxFlowAlg2(ui *ids, ui R_end, ui *pstart, ui *pend, ui *edges, int s);
    bool verify_Sbundle_by_maxFlowAlg2(ui *ids, ui S_end, ui R_end, ui *pstart, ui *pend, ui *edges, int s, ui u, ui lb);
    bool verify_Sbundle_by_maxFlowAlg3(std::vector<ui> ids, ui R_end, ui *pstart, ui *pend, ui *edges, int s);
    bool is_connected_graph(ui *ids, ui R_end, char *matrix, long long row);
    bool is_connected_graph(std::vector<ui> ids, ui R_end, char *matrix, long long row);
    bool verify_kVC_by_maxFlowAlg(ui *ids, ui R_end, ui *pstart, ui *pend, ui *edges, int s);
};
#endif