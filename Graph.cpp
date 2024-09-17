#include "Graph.h"
#include "SBundle_BB_matrix.h"
#include "SBundle_BB.h"
#include "CTPrune.h"
#include <fstream>

using namespace std;
ui pre_n;
ui pre_m;

Graph::Graph(const char *_dir, const int _K)
{
	dir = string(_dir);
	K = _K;
	S = 0;
	n = m = 0;

	pstart = nullptr;
	pend = pend_buf = nullptr;
	edges = nullptr;
	edgelist_pointer = nullptr;

	sbundle.clear();
}

Graph::~Graph()
{
	if (pstart != nullptr)
	{
		delete[] pstart;
		pstart = nullptr;
	}
	if (pend != nullptr)
	{
		delete[] pend;
		pend = nullptr;
	}
	if (pend_buf != nullptr)
	{
		delete[] pend_buf;
		pend_buf = nullptr;
	}
	if (edges != nullptr)
	{
		delete[] edges;
		edges = nullptr;
	}
	if (edgelist_pointer != nullptr)
	{
		delete[] edgelist_pointer;
		edgelist_pointer = nullptr;
	}
}


void Graph::readBinFile()
{
	ifstream in(dir);
	if (!in.is_open())
	{
		printf("Failed to open %s \n", dir.c_str());
		fflush(stdout);
		exit(1);
	}
	string suffix = "bin";
	if (suffix == "bin")
	{
		FILE *in = fopen(dir.c_str(), "rb");
		if (in == nullptr)
		{
			printf("Failed to open %s \n", dir.c_str());
			exit(1);
		}
		ui size_int;
		fread(&size_int, sizeof(ui), 1, in);
		if (size_int != sizeof(ui))
		{
			printf("sizeof int is different: graph_file(%u), machine(%u)\n", size_int, (int)sizeof(ui));
			exit(1);
		}
		fread(&n, sizeof(ui), 1, in);
		fread(&m, sizeof(ui), 1, in);
		cout << "File: " << get_file_name_without_suffix(dir) << " n= " << n << " m= " << m / 2 << " s = " << K << endl;
		ui *degree = new ui[n];
		if (pstart == nullptr)
			pstart = new ept[n + 1];
		if (edges == nullptr)
			edges = new ui[m];
		fread(degree, sizeof(ui), n, in);
		fread(edges, sizeof(ui), m, in);
		pstart[0] = 0;
		for (ui i = 1; i <= n; i++)
			pstart[i] = pstart[i - 1] + degree[i - 1];
		delete[] degree;
	}

	printf("Graph init ok\n");
	fflush(stdout);
	in.close();
}

std::string Graph::get_file_name_without_suffix(const std::string &file_path)
{
	size_t last_slash = file_path.find_last_of("/\\");
	size_t last_dot = file_path.find_last_of(".");
	return file_path.substr(last_slash + 1, last_dot - last_slash - 1);
}

int Graph::get_min_kVC()
{
	Timer t;
	long long branch_node = 0;

	assert(K > 0);
	if (K >= n)
	{
		printf("Note! There is no K-VC (K >= n) in the graph! \n");
		return 0;
	}

	pre_n = n;
	pre_m = m;
	ui *peel_sequence = new ui[n];
	ui *core = new ui[n];
	ui *degree = new ui[n];
	char *vis = new char[n];
	ListLinearHeap *heap = new ListLinearHeap(n, n - 1);
	ui UB = degen(n, peel_sequence, core, pstart, edges, degree, vis, heap, true);
	ui LB = 0;
	ui old_size = sbundle.size();
	ui *out_mapping = new ui[n];
	ui *rid = new ui[n];

	kVC_core_shrink_graph(n, m, peel_sequence, core, out_mapping, nullptr, rid, pstart, edges, true);

	if (n <= K)
	{
		printf("No K-VC!!! The whole graph is reduced to empty by core. Total Time: %s (microseconds)\n", Utility::integer_to_string(t.elapsed()).c_str());
		return -1;
	}

	vector<ui> heu_kVC;
	bool is_found = false;
	ui upper_bound = n;

	get_heuristic_k_VC(heu_kVC, is_found, upper_bound); // computes a heuritic k-VC and use its size as an upper bound of the minimum k-VC

	if (is_found) // We terminate the algorithm if the found heuristic solution is indeed the minimum k-VC, i.e., its size is k+1
	{
		assert(heu_kVC.size() > K);
		sbundle.clear();
		for (ui i = 0; i < K + 1; i++)
			sbundle.push_back(heu_kVC[i]);
		printf("The heuristic solution is indeed the minimum k-VC!!! \n");
		printf("Next print the vertex set of the obtained minimum k-VC: \n");
		std::sort(sbundle.begin(), sbundle.end()); 
		for (ui num : sbundle)
			printf("%u, ", num);
		printf("\n");
		printf("\tMinimum %d-VC Size: %lu, Total Time: %s (microseconds)\n", K, sbundle.size(), Utility::integer_to_string(t.elapsed()).c_str());
		return heu_kVC.size();
	}
	else // Otherwise, if the heuristic is not a minimum k-VC, we enter the precise search, which corresponds to Lines 3-12 in Algorithm 2 in our paper 
	{
		assert(n > K);
		Timer tt;
		SBundle_BB *sbundle_solver = new SBundle_BB();
		sbundle_solver->allocateMemory(n);

		pre_n = n, pre_m = m;
		std::vector<ui> peel_sequence_bk(peel_sequence, peel_sequence + n); 
		std::vector<ui> core_bk(core, core + n);
		std::vector<ui> degree_bk(degree, degree + n);
		std::vector<ui> out_mapping_bk(out_mapping, out_mapping + n);
		std::vector<ui> rid_bk(rid, rid + n);
		std::vector<ui> pstart_bk(pstart, pstart + n + 1);
		std::vector<ui> edges_bk(edges, edges + m);
		assert(is_found == false);

		assert(n > K && S == 0);
		for (ui ii = 0; ii < pre_n - K; ++ii) // the maximum number of the procedure solve-FSVC
		{
			assert(n == pre_n && m == pre_m); 
			S++;							  
			LB = K + S; 

			if (LB == upper_bound) // means the heuristic solution is the minimum k-VC
			{
				sbundle = heu_kVC;
				is_found = true;
				break; 
			}

			sbundle.clear(); 
			vector<ui> ids;
			vector<pair<ui, ui>> vp;

			ui *peel_sequence_rid = new ui[n];
			for (ui i = 0; i < n; i++)
				peel_sequence_rid[peel_sequence[i]] = i; 

			memset(vis, 0, sizeof(char) * n);
			SBUNDLE_BB_matrix *sbundle_solver_m = new SBUNDLE_BB_matrix(); 
			sbundle_solver_m->allocateMemory(n);

			if (pend == nullptr)
				pend = new ept[n + 1];

			reorganize_adjacency_lists(n, peel_sequence, rid, pstart, pend, edges);
			assert(sbundle.size() == 0);
			for (ui i = n; i > 0 && sbundle.size() < LB; i--) 
			{
				ui u = peel_sequence[i - 1];											 
				if (pend[u] - pstart[u] + S <= sbundle.size() || n - i < sbundle.size()) 
					continue;

				fflush(stdout);
				
				// extract u and its 2-hop neighbors
				if (sbundle.size() >= 2 * S - 1) 
					extract_subgraph_with_prune(u, sbundle.size() + 1 - S, sbundle.size() + 1 - 2 * S, sbundle.size() + 3 - 2 * S, peel_sequence_rid, degree, ids, rid, vp, vis, pstart, pend, edges);
				else
					extract_subgraph_wo_prune(u, peel_sequence_rid, ids, rid, vp, vis, pstart, pend, edges);

				if (ids.empty() || ids.size() <= sbundle.size())
					continue;

				ui t_old_size = sbundle.size();
				sbundle_solver_m->load_graph(ids.size(), vp);
				sbundle_solver_m->sBundle(S, sbundle.size() + 1, sbundle, true, vp.size(), branch_node);

				if (sbundle.size() >= LB)
				{
					is_found = true;
					for (ui j = 0; j < sbundle.size(); j++)
						sbundle[j] = ids[sbundle[j]];
					for (ui j = 0; j < sbundle.size(); j++)
						sbundle[j] = out_mapping[sbundle[j]];
					break;
				}
			}

			delete sbundle_solver_m;
			if (peel_sequence_rid != nullptr)
			{
				delete[] peel_sequence_rid;
				peel_sequence_rid = nullptr;
			}

			if (is_found)
				break;

			assert(!is_found); 
			
			// if the divide-and-conquer technique fails to find a large enough solution, then we search on the entire graph
			if (n > sbundle.size() && sbundle.size() < 2 * S - 2 && 2 * S - 2 >= LB) 
			{
				if (sbundle.size() > old_size) 
				{
					old_size = sbundle.size(); 
					core_shrink_graph(n, m, peel_sequence, core, out_mapping, nullptr, rid, pstart, edges, true);
				}
				
				sbundle_solver->load_graph(n, pstart, pstart + 1, edges);	 
				sbundle_solver->sBundle(S, LB, sbundle, false, branch_node); 
				if (sbundle.size() > 2 * S - 2)
					printf("!!! WA in sBundle_exact!\n"); 
				if (sbundle.size() == LB)
				{
					is_found = true;
					for (ui i = 0; i < sbundle.size(); i++)
						sbundle[i] = out_mapping[sbundle[i]]; 
				}
			}

			if (is_found)
				break;

			n = pre_n;
			m = pre_m;
			for (ui i = 0; i < n; i++)
			{
				peel_sequence[i] = peel_sequence_bk[i];
				core[i] = core_bk[i];
				degree[i] = degree_bk[i];
				out_mapping[i] = out_mapping_bk[i];
				rid[i] = rid_bk[i];
			}
			for (ui i = 0; i < n + 1; i++)
				pstart[i] = pstart_bk[i];
			for (ui i = 0; i < m; i++)
				edges[i] = edges_bk[i];
		}

		delete sbundle_solver;
		if (is_found)
		{
			
			printf("\tMinimum %d-VC Size: %lu, Total Time: %s (microseconds)\n", K, sbundle.size(), Utility::integer_to_string(t.elapsed()).c_str());
			printf("Next print the vertex set of the obtained minimum k-VC: \n");
			std::sort(sbundle.begin(), sbundle.end()); 
			for (ui num : sbundle)
				printf("%u, ", num);
			printf("\n");
			return LB; 
		}
		else
		{
			printf("There is no K-VC at all !!! Total Time: %s (microseconds)\n", Utility::integer_to_string(t.elapsed()).c_str());
			return -1; 
		}
	}

	delete[] out_mapping;
	delete[] rid;
	delete heap;
	delete[] core;
	delete[] peel_sequence;
	delete[] vis;
	delete[] degree;
}



void Graph::get_heuristic_k_VC(std::vector<ui> &heu_kVC, bool &is_found, ui &upper_bound)
{
	Timer t;
	ui threshold = 0;
	ui *id_s = new ui[n];
	ui *degree = new ui[n];
	ui *core = new ui[n];
	char *vis = new char[n];
	memset(vis, 0, sizeof(char) * n);

	for (ui i = 0; i < n; i++)
	{
		degree[i] = pstart[i + 1] - pstart[i];
	}

	ui queue_n = 0, new_size = 0;
	for (ui i = 0; i < n; i++)
		if (degree[i] < threshold)
			id_s[queue_n++] = i;
	for (ui i = 0; i < queue_n; i++)
	{
		ui u = id_s[i];
		degree[u] = 0;
		for (ui j = pstart[u]; j < pstart[u + 1]; j++)
			if (degree[edges[j]] > 0)
			{
				if ((degree[edges[j]]--) == threshold)
					id_s[queue_n++] = edges[j];
			}
	}
	for (ui i = 0; i < n; i++)
	{
		if (degree[i] >= threshold)
			id_s[queue_n + (new_size++)] = i;
		else
		{
			vis[i] = 1;
			core[i] = 0;
		}
	}
	assert(queue_n + new_size == n);
	ListLinearHeap *heap = new ListLinearHeap(n, new_size - 1);
	heap->init(new_size, new_size - 1, id_s + queue_n, degree);
	ui max_core = 0;
	ui id_clique = n;
	for (ui i = 0; i < new_size; i++)
	{
		ui u, key;
		heap->pop_min(u, key);
		if (key > max_core)
			max_core = key;
		core[u] = max_core;
		id_s[queue_n + i] = u;
		if (key + i + 1 == new_size)
		{
			if (id_clique == n)
				id_clique = i;
			ui x_size = i + 1;
			heap->get_ids(id_s + queue_n, x_size);
			assert(x_size == new_size);
			for (ui j = i; j < new_size; j++)
			{
				core[id_s[queue_n + j]] = max_core;
				heu_kVC.pb(id_s[queue_n + j]);
			}
			break;
		}
		vis[u] = 1;

		for (ui j = pstart[u]; j < pstart[u + 1]; j++)
			if (vis[edges[j]] == 0)
			{
				heap->decrement(edges[j], 1);
			}
	}

	if (heu_kVC.size() < K + 1)
	{
		MaximumFlow *vc = new MaximumFlow();
		vc->prepare(n, m);
		int add_vertex_idx = id_clique - 1;
		while (heu_kVC.size() < K + 1 && add_vertex_idx > queue_n)
		{
			heu_kVC.push_back(id_s[add_vertex_idx]);
			add_vertex_idx--;
		}
		assert(heu_kVC.size() >= K + 1 || add_vertex_idx <= queue_n);
		if (heu_kVC.size() >= K + 1)
		{
			assert(heu_kVC.size() == K + 1);
			bool test_kVC = vc->verify_kVC_by_maxFlowAlg(heu_kVC.data(), heu_kVC.size(), pstart, pstart + 1, edges, K);
			while (!test_kVC)
			{
				heu_kVC.push_back(id_s[add_vertex_idx]);
				add_vertex_idx--;
				if (add_vertex_idx <= queue_n)
					break;
				test_kVC = vc->verify_kVC_by_maxFlowAlg(heu_kVC.data(), heu_kVC.size(), pstart, pstart + 1, edges, K);
			}

			if (test_kVC)
			{

				upper_bound = heu_kVC.size();
				assert(upper_bound > K);
			}
			else
			{
				upper_bound = n;
			}
		}
		delete vc;
	}
	else
	{
		is_found = true;
		upper_bound = K + 1;
	}

	delete heap;
	delete[] id_s;
	delete[] core;
	delete[] degree;
	delete[] vis;
	printf("After the heuristic operation, we find the upper bound of the minium kVC is %d, and the time cost is %s\n", upper_bound, Utility::integer_to_string(t.elapsed()).c_str());
}

void Graph::reorganize_adjacency_lists(ui n, ui *peel_sequence, ui *rid, ui *pstart, ui *pend, ui *edges)
{
	for (ui i = 0; i < n; i++)
		rid[peel_sequence[i]] = i;

	for (ui i = 0; i < n; i++)
	{
		ui &end = pend[i] = pstart[i];
		for (ui j = pstart[i]; j < pstart[i + 1]; j++)
			if (rid[edges[j]] > rid[i])
				edges[end++] = edges[j];
	}

	for (ui i = n; i > 0; i--)
	{
		ui u = peel_sequence[i - 1];
		for (ui j = pstart[u]; j < pend[u] && rid[edges[j]] > rid[u]; j++)
		{
			ui v = edges[j];
			edges[pend[v]++] = u;
			assert(pend[v] <= pstart[v + 1]);
		}
	}
#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(pend[i] == pstart[i + 1]);
#endif
	for (ui i = 0; i < n; i++)
	{
		ui &end = pend[i] = pstart[i];
		while (end < pstart[i + 1] && rid[edges[end]] > rid[i])
			++end;
	}
}

void Graph::extract_subgraph_with_prune(ui u, ui degree_threshold, ui triangle_threshold, ui cn_threshold, const ui *p_rid, ui *degree, vector<ui> &ids, ui *rid, vector<pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges)
{
#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(!exists[i]);
#endif

	ids.clear();
	vp.clear();
	ids.push_back(u);
	exists[u] = 1;
	for (ept i = pstart[u]; i < pend[u]; i++)
	{
		assert(p_rid[edges[i]] > p_rid[u]);
		ids.push_back(edges[i]);
		exists[edges[i]] = 2;
	}
	assert(pend[u] >= pstart[u + 1] || p_rid[edges[pend[u]]] < p_rid[u]);

	ui *Q = rid;
	ui Q_n = 0;
	for (ui i = 1; i < ids.size(); i++)
	{
		ui v = ids[i];
		degree[v] = 0;
		for (ept j = pstart[v]; j < pstart[v + 1] && p_rid[edges[j]] > p_rid[u]; j++)
		{
			if (exists[edges[j]])
				++degree[v];
		}
		if (degree[v] < triangle_threshold)
			Q[Q_n++] = v;
	}

	for (ui i = 0; i < Q_n; i++)
	{
		ui v = Q[i];
		exists[v] = 3;
		for (ept j = pstart[v]; j < pstart[v + 1] && p_rid[edges[j]] > p_rid[u]; j++)
			if (exists[edges[j]] == 2)
			{
				if (degree[edges[j]] == triangle_threshold)
					Q[Q_n++] = edges[j];
				--degree[edges[j]];
			}
	}
	assert(Q_n < ids.size());

	if (ids.size() - Q_n - 1 < degree_threshold)
	{
		for (ui i = 0; i < ids.size(); i++)
			exists[ids[i]] = 0;
		ids.clear();
		return;
	}

	ui old_size = ids.size();
	for (ui i = 1; i < old_size; i++)
		if (exists[ids[i]] == 2)
		{
			ui v = ids[i];
			for (ept j = pstart[v]; j < pstart[v + 1] && p_rid[edges[j]] > p_rid[u]; j++)
			{
				if (!exists[edges[j]])
				{
					ids.push_back(edges[j]);
					exists[edges[j]] = 1;
					degree[edges[j]] = 1;
				}
				else
					++degree[edges[j]];
			}
		}

	ui new_size = 1;
	for (ui i = 1; i < old_size; i++)
	{
		if (exists[ids[i]] == 3)
			exists[ids[i]] = 0;
		else
			ids[new_size++] = ids[i];
	}
	assert(new_size + Q_n == old_size);

	for (ui i = old_size; i < ids.size(); i++)
	{
		if (degree[ids[i]] < cn_threshold)
			exists[ids[i]] = 0;
		else
			ids[new_size++] = ids[i];
	}
	ids.resize(new_size);

	for (ui i = 0; i < ids.size(); i++)
		rid[ids[i]] = i;

	for (ui i = 0; i < ids.size(); i++)
	{
		ui v = ids[i];
		for (ept j = pstart[v]; j < pend[v]; j++)
			if (exists[edges[j]])
			{
				assert(rid[v] < ids.size() && rid[edges[j]] < ids.size());
				vp.push_back(make_pair(rid[v], rid[edges[j]]));
			}
	}

	for (ui i = 0; i < ids.size(); i++)
		exists[ids[i]] = 0;
}

void Graph::extract_subgraph_wo_prune(ui u, const ui *p_rid, vector<ui> &ids, ui *rid, vector<pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges)
{
	ids.clear();
	vp.clear();
	ids.push_back(u);
	exists[u] = 1;
	rid[u] = 0;

	for (ept i = pstart[u]; i < pend[u]; i++)
	{
		assert(p_rid[edges[i]] > p_rid[u]);
		ids.push_back(edges[i]);
		exists[edges[i]] = 1;
		rid[edges[i]] = ids.size() - 1;
	}
	assert(pend[u] >= pstart[u + 1] || p_rid[edges[pend[u]]] < p_rid[u]);
	ui old_size = ids.size();
	for (ui i = 1; i < old_size; i++)
	{
		ui v = ids[i];
		for (ept j = pstart[v]; j < pstart[v + 1] && p_rid[edges[j]] > p_rid[u]; j++)
		{
			ui w = edges[j];
			if (exists[w])
				continue;
			ids.push_back(w);
			exists[w] = 1;
			rid[w] = ids.size() - 1;
		}
	}
	for (ui i = 0; i < ids.size(); i++)
	{
		ui v = ids[i];
		for (ept j = pstart[v]; j < pend[v]; j++)
			if (exists[edges[j]])
				vp.push_back(make_pair(rid[v], rid[edges[j]]));
	}
	for (ui i = 0; i < ids.size(); i++)
		exists[ids[i]] = 0;
}

void Graph::write_subgraph(ui n, const vector<pair<int, int>> &edge_list)
{
	FILE *fout = Utility::open_file("edges.txt", "w");

	fprintf(fout, "%u %lu\n", n, edge_list.size());
	for (ui i = 0; i < edge_list.size(); i++)
		fprintf(fout, "%d %d\n", edge_list[i].first, edge_list[i].second);

	fclose(fout);
}

void Graph::extract_subgraph_and_prune(ui u, ui *ids, ui &ids_n, ui *rid, vector<pair<ui, ui>> &vp, ui *Q, ui *degree, char *exists, ept *pend, char *deleted, ui *edgelist_pointer)
{
	vp.clear();
	ids_n = 0;
	ids[ids_n++] = u;
	exists[u] = 1;
	ui u_n = pstart[u];
	for (ept i = pstart[u]; i < pend[u]; i++)
		if (!deleted[edgelist_pointer[i]])
		{
			edges[u_n] = edges[i];
			edgelist_pointer[u_n++] = edgelist_pointer[i];
			ui v = edges[i];
			ids[ids_n++] = v;
			exists[v] = 2;
		}
	pend[u] = u_n;

	ui Q_n = 0;
	for (ui i = 1; i < ids_n; i++)
	{
		u = ids[i];
		u_n = pstart[u];
		degree[u] = 0;
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (!deleted[edgelist_pointer[j]])
			{
				edges[u_n] = edges[j];
				edgelist_pointer[u_n++] = edgelist_pointer[j];
				if (exists[edges[j]] == 2)
					++degree[u];
			}
		pend[u] = u_n;
		if (degree[u] + 2 * K <= sbundle.size())
			Q[Q_n++] = u;
	}
	for (ui i = 0; i < Q_n; i++)
	{
		u = Q[i];
		exists[u] = 10;
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (exists[edges[j]] == 2)
			{
				if ((degree[edges[j]]--) + 2 * K == sbundle.size() + 1)
				{
					assert(Q_n < m / 2);
					Q[Q_n++] = edges[j];
				}
			}
	}
	assert(Q_n <= ids_n);
	if (ids_n - 1 - Q_n + K <= sbundle.size())
	{
		for (ui i = 0; i < ids_n; i++)
			exists[ids[i]] = 0;
		ids_n = 0;
		return;
	}

	ui nr_size = ids_n;
	for (ui i = 1; i < nr_size; i++)
		if (exists[ids[i]] == 2)
		{
			u = ids[i];
			for (ept j = pstart[u]; j < pend[u]; j++)
			{
				if (!exists[edges[j]])
				{
					ids[ids_n++] = edges[j];
					exists[edges[j]] = 3;
					degree[edges[j]] = 1;
				}
				else if (exists[edges[j]] == 3)
					++degree[edges[j]];
			}
		}

#ifndef NDEBUG

#endif

	ui new_size = 1;
	for (ui i = 1; i < nr_size; i++)
	{
		if (exists[ids[i]] == 10)
			exists[ids[i]] = 0;
		else
			ids[new_size++] = ids[i];
	}
#ifndef NDEBUG
	if (new_size + Q_n != nr_size)
	{
		printf("new_size: %u, Q_n: %u, nr_size: %u\n", new_size, Q_n, nr_size);
		printf("New list: ");
		for (ui i = 0; i < new_size; i++)
			printf(" %u", ids[i]);
		printf("\n");
		printf("Pruned list: ");
		for (ui i = 0; i < Q_n; i++)
			printf(" %u", Q[i]);
		printf("\n");
	}
#endif
	assert(new_size + Q_n == nr_size);
	ui old_nr_size = nr_size;
	nr_size = new_size;
	for (ui i = old_nr_size; i < ids_n; i++)
	{
		if (degree[ids[i]] + 2 * K <= sbundle.size() + 2)
			exists[ids[i]] = 0;
		else
			ids[new_size++] = ids[i];
	}
	ids_n = new_size;
#ifndef NDEBUG
	assert(exists[ids[0]] == 1);
	for (ui i = 1; i < nr_size; i++)
		assert(exists[ids[i]] == 2);
	for (ui i = nr_size; i < ids_n; i++)
		assert(exists[ids[i]] == 3);
#endif

	for (ui i = 0; i < ids_n; i++)
	{
		assert(exists[ids[i]]);
		rid[ids[i]] = i;
	}

	for (ui i = 0; i < nr_size; i++)
	{
		u = ids[i];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (exists[edges[j]] && edges[j] > u)
			{
				assert(!deleted[edgelist_pointer[j]]);
				vp.push_back(make_pair(rid[u], rid[edges[j]]));
			}
	}
	for (ui i = nr_size; i < ids_n; i++)
	{
		u = ids[i];
		u_n = pstart[u];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (!deleted[edgelist_pointer[j]])
			{
				edges[u_n] = edges[j];
				edgelist_pointer[u_n++] = edgelist_pointer[j];
				if (edges[j] > u && exists[edges[j]])
					vp.push_back(make_pair(rid[u], rid[edges[j]]));
			}
		pend[u] = u_n;
	}
	for (ui i = 0; i < ids_n; i++)
		exists[ids[i]] = 0;
#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(exists[i] == 0);
#endif
}

ui Graph::degen(ui n, ui *peel_sequence, ui *core, ept *pstart, ui *edges, ui *degree, char *vis, ListLinearHeap *heap, bool output)
{
	Timer t;

	ui threshold = 0;
	for (ui i = 0; i < n; i++)
		degree[i] = pstart[i + 1] - pstart[i];

	ui queue_n = 0, new_size = 0;
	for (ui i = 0; i < n; i++)
		if (degree[i] < threshold)
			peel_sequence[queue_n++] = i;
	for (ui i = 0; i < queue_n; i++)
	{
		ui u = peel_sequence[i];
		degree[u] = 0;
		for (ept j = pstart[u]; j < pstart[u + 1]; j++)
			if (degree[edges[j]] > 0)
			{
				if ((degree[edges[j]]--) == threshold)
					peel_sequence[queue_n++] = edges[j];
			}
	}

	ui UB = n;
	if (queue_n == n)
		UB = sbundle.size();

	memset(vis, 0, sizeof(char) * n);
	for (ui i = 0; i < n; i++)
	{
		if (degree[i] >= threshold)
			peel_sequence[queue_n + (new_size++)] = i;
		else
		{
			vis[i] = 1;
			core[i] = 0;
		}
	}
	assert(queue_n + new_size == n);

	if (new_size != 0)
	{
		heap->init(new_size, new_size - 1, peel_sequence + queue_n, degree);
		ui max_core = 0;
		ui idx = n;
		ui idx_splex = n;
		UB = 0;
		for (ui i = 0; i < new_size; i++)
		{
			ui u, key;
			heap->pop_min(u, key);
			if (key > max_core)
				max_core = key;
			core[u] = max_core;
			peel_sequence[queue_n + i] = u;

			ui t_UB = core[u] + K;
			if (new_size - i < t_UB)
				t_UB = new_size - i;
			if (t_UB > UB)
				UB = t_UB;
			vis[u] = 1;

			for (ept j = pstart[u]; j < pstart[u + 1]; j++)
				if (vis[edges[j]] == 0)
				{
					heap->decrement(edges[j], 1);
				}
		}
	}
	return UB;
}

void Graph::core_shrink_graph(ui &n, ept &m, ui *peel_sequence, ui *core, ui *out_mapping, ui *in_mapping, ui *rid, ept *&pstart, ui *&edges, bool output)
{

	ui cnt = 0;
	for (ui i = 0; i < n; i++)
		if (core[i] + S > sbundle.size())
		{
			rid[i] = cnt;
			if (in_mapping == nullptr)
				out_mapping[cnt] = i;
			else
				out_mapping[cnt] = in_mapping[i];
			++cnt;
		}

	if (cnt != n)
	{
		cnt = 0;
		ept pos = 0;
		for (ui i = 0; i < n; i++)
			if (core[i] + S > sbundle.size())
			{
				ept t_start = pstart[i];
				pstart[cnt] = pos;
				for (ept j = t_start; j < pstart[i + 1]; j++)
					if (core[edges[j]] + S > sbundle.size())
					{
						edges[pos++] = rid[edges[j]];
					}
				++cnt;
			}
		pstart[cnt] = pos;

		assert(core[peel_sequence[n - cnt - 1]] == 0 || core[peel_sequence[n - cnt - 1]] + S <= sbundle.size());
		assert(cnt == 0 || core[peel_sequence[n - cnt]] + S > sbundle.size());

		for (ui i = 0; i < cnt; i++)
		{
			peel_sequence[i] = rid[peel_sequence[n - cnt + i]];
			core[i] = core[out_mapping[i]];
		}

		if (pos > 0 && pos < m / 2)
		{

			ept *pstart_new = new ept[cnt + 1];
			ui *edges_new = new ui[pos];
			memcpy(pstart_new, pstart, sizeof(ept) * (cnt + 1));
			memcpy(edges_new, edges, sizeof(ui) * pos);

			delete[] pstart;
			pstart = pstart_new;
			delete[] edges;
			edges = edges_new;
		}
		n = cnt;
		m = pos;
	}

	if (output)
		printf("*** After core shrink: n = %s, m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m / 2).c_str());
}

void Graph::kVC_core_shrink_graph(ui &n, ept &m, ui *peel_sequence, ui *core, ui *out_mapping, ui *in_mapping, ui *rid, ept *&pstart, ui *&edges, bool output)
{

	ui cnt = 0;
	for (ui i = 0; i < n; i++)
		if (core[i] >= K)
		{
			rid[i] = cnt;
			if (in_mapping == nullptr)
				out_mapping[cnt] = i;
			else
				out_mapping[cnt] = in_mapping[i];
			++cnt;
		}

	if (cnt != n)
	{
		cnt = 0;
		ept pos = 0;
		for (ui i = 0; i < n; i++)
			if (core[i] >= K)
			{
				ept t_start = pstart[i];
				pstart[cnt] = pos;
				for (ept j = t_start; j < pstart[i + 1]; j++)
					if (core[edges[j]] >= K)
					{
						edges[pos++] = rid[edges[j]];
					}
				++cnt;
			}
		pstart[cnt] = pos;

		assert(core[peel_sequence[n - cnt - 1]] == 0 || core[peel_sequence[n - cnt - 1]] < K);
		assert(cnt == 0 || core[peel_sequence[n - cnt]] >= K);

		for (ui i = 0; i < cnt; i++)
		{
			peel_sequence[i] = rid[peel_sequence[n - cnt + i]];
			core[i] = core[out_mapping[i]];
		}

		if (pos > 0 && pos < m / 2)
		{

			ept *pstart_new = new ept[cnt + 1];
			ui *edges_new = new ui[pos];
			memcpy(pstart_new, pstart, sizeof(ept) * (cnt + 1));
			memcpy(edges_new, edges, sizeof(ui) * pos);

			delete[] pstart;
			pstart = pstart_new;
			delete[] edges;
			edges = edges_new;
		}
		n = cnt;
		m = pos;
	}

	if (output)
		printf("*** After core shrink: n = %s, m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m / 2).c_str());
}

void Graph::orient_graph(ui n, ui m, ui *peel_sequence, ept *pstart, ept *pend, ui *edges, ui *rid)
{
	for (ui i = 0; i < n; i++)
		rid[peel_sequence[i]] = i;
	for (ui i = 0; i < n; i++)
	{
		ept &end = pend[i] = pstart[i];
		for (ept j = pstart[i]; j < pstart[i + 1]; j++)
			if (rid[edges[j]] > rid[i])
				edges[end++] = edges[j];
	}

#ifndef NDEBUG
	long long sum = 0;
	for (int i = 0; i < n; i++)
		sum += pend[i] - pstart[i];
	assert(sum * 2 == m);
#endif
}

void Graph::oriented_triangle_counting(ui n, ui m, ept *pstart, ept *pend, ui *edges, ui *tri_cnt, ui *adj)
{
	memset(adj, 0, sizeof(ui) * n);
	long long cnt = 0;
	memset(tri_cnt, 0, sizeof(ui) * m);
	for (ui u = 0; u < n; u++)
	{
		for (ept j = pstart[u]; j < pend[u]; j++)
			adj[edges[j]] = j + 1;

		for (ept j = pstart[u]; j < pend[u]; j++)
		{
			ui v = edges[j];
			for (ept k = pstart[v]; k < pend[v]; k++)
				if (adj[edges[k]])
				{
					++tri_cnt[j];
					++tri_cnt[k];
					++tri_cnt[adj[edges[k]] - 1];
					++cnt;
				}
		}

		for (ept j = pstart[u]; j < pend[u]; j++)
			adj[edges[j]] = 0;
	}

#ifndef NDEBUG

#endif
}

void Graph::reorganize_oriented_graph(ui n, ui *tri_cnt, ui *edge_list, ept *pstart, ept *pend, ept *pend2, ui *edges, ui *edgelist_pointer, ui *buf)
{
	for (ui i = 0; i < n; i++)
		pend2[i] = pend[i];
	ept pos = 0;
	for (ui i = 0; i < n; i++)
	{
		for (ept j = pstart[i]; j < pend[i]; j++)
		{
			tri_cnt[pos >> 1] = edgelist_pointer[j];
			edge_list[pos++] = i;
			edge_list[pos++] = edges[j];

			ept &k = pend2[edges[j]];
			edgelist_pointer[k] = edgelist_pointer[j] = (pos >> 1) - 1;
			edges[k++] = i;
		}
	}

#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(pend2[i] == pstart[i + 1]);
#endif

	for (ui i = 0; i < n; i++)
	{
		pend2[i] = pend[i];
		pend[i] = pstart[i];
	}
	for (ui i = 0; i < n; i++)
	{
		for (ept j = pend2[i]; j < pstart[i + 1]; j++)
		{
			ept &k = pend[edges[j]];
			edgelist_pointer[k] = edgelist_pointer[j];
			edges[k++] = i;
		}
	}

	ept *ids = pend2;
	for (ui i = 0; i < n; i++)
	{
		if (pend[i] == pstart[i] || pend[i] == pstart[i + 1])
			continue;
		ept j = pstart[i], k = pend[i], pos = 0;
		while (j < pend[i] && k < pstart[i + 1])
		{
			if (edges[j] < edges[k])
			{
				ids[pos] = edges[j];
				buf[pos++] = edgelist_pointer[j++];
			}
			else
			{
				ids[pos] = edges[k];
				buf[pos++] = edgelist_pointer[k++];
			}
		}
		while (j < pend[i])
		{
			ids[pos] = edges[j];
			buf[pos++] = edgelist_pointer[j++];
		}
		while (k < pstart[i + 1])
		{
			ids[pos] = edges[k];
			buf[pos++] = edgelist_pointer[k++];
		}
		for (ept j = 0; j < pos; j++)
		{
			edges[pstart[i] + j] = ids[j];
			edgelist_pointer[pstart[i] + j] = buf[j];
		}
	}
}
