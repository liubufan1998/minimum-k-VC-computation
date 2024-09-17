#ifndef _SBUNDLE_BB_
#define _SBUNDLE_BB_

#include "Utility.h"
#include "Timer.h"
#include "sbundle_tool.h"

class SBundle_BB
{
private:
	ui n;
	ui temp_n;
	ept *pstart;
	ept *pend;
	ui *edges;
	ept edges_cap;

	ui *degree;
	std::vector<ui> degree2;
	ui *degree_in_S;

	ui s;
	ui *best_solution;
	ui best_solution_size;
	ui _UB_;

	ui *neighbors;
	ui *nonneighbors;

	ui *S2;
	ui *SR;
	ui *SR_rid;
	ui *SR_remap;
	ui *SR_remap_rid;
	std::queue<ui> Qv;
	ui *level_id;
	ui *t_level_id;

	ui *buf;
	ui *buf1;
	ui *buf2;
	char *vis;

	std::vector<std::pair<ui, ui>> vp;
	char *matrix;
	std::vector<ui> temp_deleted;
	MaximumFlow *vc;
	long long branch_count;
	std::vector<ui> temp_test;

public:
	SBundle_BB()
	{
		n = 0;
		edges_cap = 0;
		pstart = NULL;
		pend = NULL;
		edges = NULL;
		degree = degree_in_S = NULL;

		best_solution = NULL;
		s = best_solution_size = _UB_ = 0;

		S2 = SR = SR_rid = NULL;
		level_id = NULL;

		neighbors = nonneighbors = NULL;
		buf = buf1 = buf2 = NULL;
		vis = NULL;
	}

	~SBundle_BB()
	{
		if (pstart != NULL)
		{
			delete[] pstart;
			pstart = NULL;
		}
		if (pend != NULL)
		{
			delete[] pend;
			pend = NULL;
		}
		if (edges != NULL)
		{
			delete[] edges;
			edges = NULL;
		}
		if (degree != NULL)
		{
			delete[] degree;
			degree = NULL;
		}
		if (degree_in_S != NULL)
		{
			delete[] degree_in_S;
			degree_in_S = NULL;
		}
		if (best_solution != NULL)
		{
			delete[] best_solution;
			best_solution = NULL;
		}
		if (S2 != NULL)
		{
			delete[] S2;
			S2 = NULL;
		}
		if (SR != NULL)
		{
			delete[] SR;
			SR = NULL;
		}
		if (SR_rid != NULL)
		{
			delete[] SR_rid;
			SR_rid = NULL;
		}
		if (level_id != NULL)
		{
			delete[] level_id;
			level_id = NULL;
		}
		if (neighbors != NULL)
		{
			delete[] neighbors;
			neighbors = NULL;
		}
		if (nonneighbors != NULL)
		{
			delete[] nonneighbors;
			nonneighbors = NULL;
		}
		if (buf != NULL)
		{
			delete[] buf;
			buf = NULL;
		}
		if (buf1 != NULL)
		{
			delete[] buf1;
			buf1 = NULL;
		}
		if (buf2 != NULL)
		{
			delete[] buf2;
			buf2 = NULL;
		}
		if (vis != NULL)
		{
			delete[] vis;
			vis = NULL;
		}
	}

	void allocateMemory(ui n)
	{
		if (n <= 0)
			return;
		pstart = new ept[n + 1];
		pend = new ept[n + 1];
		edges = new ui[1];
		edges_cap = 1;

		degree = new ui[n + 1];
		degree_in_S = new ui[n + 1];
		best_solution = new ui[n + 1];
		S2 = new ui[n + 1];
		SR = new ui[n + 1];
		SR_rid = new ui[n + 1];
		neighbors = new ui[n + 1];
		nonneighbors = new ui[n + 1];
		level_id = new ui[n + 1];

		buf = new ui[n + 1];
		buf1 = new ui[n + 1];
		buf2 = new ui[n + 1];
		vis = new char[n + 1];
	}

	void load_graph(ui n_, ui *_pstart, ui *_pend, ui *_edges)
	{
		ept m = 0;
		n = n_;
		for (ui i = 0; i < n; i++)
			m += _pend[i] - _pstart[i];
		if (m > edges_cap)
		{
			do
			{
				edges_cap *= 2;
			} while (m > edges_cap);
			delete[] edges;
			edges = new ui[edges_cap];
		}

		m = 0;
		for (ui i = 0; i < n; i++)
		{
			pstart[i] = m;
			for (ept j = _pstart[i]; j < _pend[i]; j++)
			{
				assert(_edges[j] < n);
				edges[m++] = _edges[j];
			}
		}
		pstart[n] = m;

		printf("load graph of size n=%u, m=%u (undirected), density=%.5lf\n", n, m / 2, double(m) / n / (n - 1));
	}

	void sBundle(ui S_, ui UB_, std::vector<ui> &sbundle, bool must_include_0, long long &branch_node)
	{
		s = S_;
		_UB_ = UB_;
		branch_count = 0;

		best_solution_size = sbundle.size();
		ui R_end;
		initialization(R_end);
		vc = new MaximumFlow();
		vc->prepare(n, edges_cap);

		if (R_end && best_solution_size < _UB_)
			splex_BB_search(0, R_end, 1, must_include_0);

		branch_node += branch_count;
		if (vc != nullptr)
		{
			delete vc;
			vc = nullptr;
		}
		if (best_solution_size > sbundle.size())
		{
			sbundle.clear();
			for (ui i = 0; i < best_solution_size; i++)
				sbundle.push_back(best_solution[i]);
		}
	}

	bool canadd(ui u, ui S_end, bool precise)
	{
		char *t_matrix = matrix + u * temp_n;

		if (precise)
		{
			if (SR_remap_rid[u] != S_end)
			{
				swap_pos2(S_end, SR_remap_rid[u]);
			}
			return vc->verify_SBundle_by_MaxFlowAlg0(SR_remap, S_end + 1, matrix, temp_n, s, false);
		}

		return true;
	}

	void delfrC(ui u, ui &R_end, ui level)
	{
		assert(t_level_id[u] == temp_n);
		t_level_id[u] = level;
		--R_end;
		swap_pos2(R_end, SR_remap_rid[u]);

		char *t_matrix = matrix + u * temp_n;
		for (ui i = 0; i < R_end; i++)
			if (t_matrix[SR_remap[i]])
			{
				assert(degree2[SR_remap[i]] > 0);
				--degree2[SR_remap[i]];
			}
	}

	void addtoC(ui u, ui &R_end, ui level)
	{
		assert(t_level_id[u] != temp_n && t_level_id[u] == level);
		t_level_id[u] = temp_n;
		SR_remap_rid[u] = R_end;
		SR_remap[R_end] = u;
		++R_end;

		char *t_matrix = matrix + u * temp_n;
		degree2[u] = 0;
		for (ui i = 0; i < R_end; i++)
		{
			if (t_matrix[SR_remap[i]])
			{
				++degree2[SR_remap[i]];
				degree2[u]++;
			}
		}
	}

	void addtoS(ui u, ui &S_end, ui R_end)
	{
		swap_pos2(S_end, SR_remap_rid[u]);
		++S_end;
	}

	void delfrS(ui u, ui &S_end, ui R_end)
	{
		--S_end;
		assert(S_end == SR_remap_rid[u]);
	}

	void initialization(ui &R_end)
	{
		memset(vis, 0, sizeof(char) * (n + 1));
		memset(degree_in_S, 0, sizeof(ui) * (n + 1));
		for (ui i = 0; i < n; i++)
			level_id[i] = n;
		for (ui i = 0; i < n; i++)
			SR[i] = SR_rid[i] = i;
		for (ui i = 0; i < n; i++)
		{
			pend[i] = pstart[i + 1];
		}

		while (!Qv.empty())
			Qv.pop();

		for (ui i = 0; i < n; i++)
		{
			degree[i] = pstart[i + 1] - pstart[i];
			if (degree[i] + s <= best_solution_size)
			{
				level_id[i] = 0;
				Qv.push(i);
			}
		}
		R_end = n;
		if (!remove_vertices_with_prune(0, R_end, 0))
			R_end = 0;
	}

	void reorganize_edges(ui S_end, ui R_end, ui level)
	{
		assert(level > 0);
		for (ui i = 0; i < R_end; i++)
		{
			ui u = SR[i];
			assert(level_id[u] > level);
			ui non_neighbors_n = 0, end = pstart[u];
			for (ui j = pstart[u]; j < pend[u] && level_id[edges[j]] >= level - 1; j++)
			{
				assert(level_id[edges[j]] == level - 1 || level_id[edges[j]] == n);
				if (level_id[edges[j]] >= level)
				{
					edges[end++] = edges[j];
				}
				else
					nonneighbors[non_neighbors_n++] = edges[j];
			}

			assert(degree[u] == end - pstart[u]);
			for (ui j = 0; j < non_neighbors_n; j++)
				edges[end++] = nonneighbors[j];
			assert((end < pend[u] && level_id[edges[end]] < level - 1) || end == pend[u]);
#ifndef NDEBUG
			for (ui j = end; j < pend[u]; j++)
			{
				if (level_id[edges[j]] >= level)
					printf("removed_level[edges[j]]: %u, level: %u\n", level_id[edges[j]], level);
				assert(level_id[edges[j]] < level);
			}
#endif
		}
	}

	void store_solution(ui size)
	{
		if (size <= best_solution_size)
		{
			printf("!!! the solution to store is no larger than the current best solution!");
			return;
		}
		best_solution_size = size;
		for (ui i = 0; i < best_solution_size; i++)
			best_solution[i] = SR[i];
	}

	bool is_kplex(ui R_end)
	{
		for (ui i = 0; i < R_end; i++)
			if (degree[SR[i]] + s < R_end)
				return false;
		return true;
	}

	void splex_BB_search(ui S_end, ui R_end, ui level, bool choose_zero)
	{

		if (S_end > best_solution_size && vc->verify_Sbundle_by_maxFlowAlg2(SR, S_end, pstart, pend, edges, s))
			store_solution(S_end);
		if (R_end > best_solution_size && is_kplex(R_end) && vc->verify_Sbundle_by_maxFlowAlg2(SR, R_end, pstart, pend, edges, s))
			store_solution(R_end);
		if (R_end <= best_solution_size + 1 || best_solution_size >= _UB_)
			return;

		branch_count++;

#ifndef NDEBUG
		for (ui i = 0; i < S_end; i++)
		{
			assert(degree[SR[i]] + s > best_solution_size);
			assert(degree_in_S[SR[i]] + s >= S_end);
		}
#endif

		ui old_S_end = S_end, old_R_end = R_end;
		assert(Qv.empty());

		if (level > 1)
			reorganize_edges(S_end, R_end, level);
#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui u = SR[i], cnt = 0;
			for (ui j = pstart[u]; j < pend[u] && level_id[edges[j]] >= level; j++)
				if (SR_rid[edges[j]] < R_end)
				{
					++cnt;
				}
			assert(degree[u] == cnt);
			assert(level_id[u] > level);
		}
#endif
		assert(choose_zero == false);

		if (R_end > best_solution_size && is_kplex(R_end) && vc->verify_Sbundle_by_maxFlowAlg2(SR, R_end, pstart, pend, edges, s))
			store_solution(R_end);
		if (R_end <= best_solution_size + 1 || best_solution_size >= _UB_)
		{
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
			return;
		}

		ui u = n;
		ui min_deg = n;
		for (ui i = 0; i < R_end; i++)
		{
			if (degree[SR[i]] < min_deg)
			{
				u = SR[i];
				min_deg = degree[SR[i]];
			}
		}
		assert(u != n && SR_rid[u] < R_end);

		if (min_deg >= R_end - s)
		{

			ui pre_S_end = S_end, pre_R_end = R_end;
			sbundle_binary_BB(S_end, R_end, level + 1);
			assert(S_end == pre_S_end && R_end == pre_R_end);
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
			return;
		}

		assert(min_deg < R_end - s);

		if (SR_rid[u] < S_end)
		{
			u = choose_branch_vertex_based_on_non_neighbors(S_end, R_end, u);
		}
		assert(SR_rid[u] >= S_end && SR_rid[u] < R_end);
		assert(degree[u] + s > best_solution_size && degree[u] + s > S_end);

		bool can_add_u = false;
		if (S_end + 1 <= s)
		{
			can_add_u = true;
		}
		else
		{
			swap_pos(S_end, SR_rid[u]);
			can_add_u = vc->verify_Sbundle_by_maxFlowAlg2(SR, S_end + 1, pstart, pend, edges, s);
		}

		if (can_add_u)
		{

			ui pre_best_solution_size = best_solution_size, t_old_S_end = S_end, t_old_R_end = R_end;
			if (move_u_to_S_with_prune(u, S_end, R_end, level))
			{
				splex_BB_search(S_end, R_end, level + 1, false);
			}
			if (best_solution_size >= _UB_)
				return;
			assert(S_end == t_old_S_end + 1 && SR[S_end - 1] == u);
			restore_SR(S_end, R_end, S_end, t_old_R_end, level);

#ifndef NDEBUG
			for (ui i = 0; i < n; i++)
				assert(!vis[i]);
#endif

			assert(Qv.empty());
			bool succeed = remove_u_from_S_with_prune(S_end, R_end, level);
			if (succeed && best_solution_size > pre_best_solution_size)
				succeed = collect_removable_vertices(S_end, R_end, level);
			if (succeed)
				succeed = remove_vertices_with_prune(S_end, R_end, level);
			if (succeed)
			{
				splex_BB_search(S_end, R_end, level + 1, false);
			}
			if (best_solution_size >= _UB_)
				return;
			assert(S_end >= old_S_end && R_end <= old_R_end);
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
		}
		else
		{
			ui pre_best_solution_size = best_solution_size;

			assert(Qv.empty());
			bool succeed = remove_u_from_C_with_prune(u, S_end, R_end, level);
			if (succeed && best_solution_size > pre_best_solution_size)
				succeed = collect_removable_vertices(S_end, R_end, level);
			if (succeed)
				succeed = remove_vertices_with_prune(S_end, R_end, level);
			if (succeed)
			{
				splex_BB_search(S_end, R_end, level + 1, false);
			}
			if (best_solution_size >= _UB_)
				return;
			assert(S_end >= old_S_end && R_end <= old_R_end);
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
		}
	}

	void sbundle_binary_BB(ui S_end, ui R_end, ui level)
	{
		branch_count++;

		if (S_end > best_solution_size && vc->verify_Sbundle_by_maxFlowAlg2(SR, S_end, pstart, pend, edges, s))
			store_solution(S_end);
		if (R_end > best_solution_size && is_kplex(R_end) && vc->verify_Sbundle_by_maxFlowAlg2(SR, R_end, pstart, pend, edges, s))
			store_solution(R_end);
		if (R_end <= best_solution_size + 1 || best_solution_size >= _UB_)
			return;

		if (level > 1)
			reorganize_edges(S_end, R_end, level);
#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui u = SR[i], cnt = 0;
			for (ui j = pstart[u]; j < pend[u] && level_id[edges[j]] >= level; j++)
				if (SR_rid[edges[j]] < R_end)
				{
					++cnt;
				}
			assert(degree[u] == cnt);
			assert(level_id[u] > level);
		}
#endif

		ui u = n;
		ui min_deg = n;

		for (ui i = S_end; i < R_end; i++)
		{
			if (degree[SR[i]] < min_deg)
			{
				u = SR[i];
				min_deg = degree[SR[i]];
			}
		}
		assert(u != n && min_deg != n);
		assert(SR_rid[u] < R_end && SR_rid[u] >= S_end);
		assert(level_id[u] == n);
		level_id[u] = level;
		--R_end;
		swap_pos(R_end, SR_rid[u]);

		for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
		{
			if (SR_rid[edges[i]] < R_end)
			{
				ui w = edges[i];
				--degree[w];
			}
		}

		sbundle_binary_BB(S_end, R_end, level + 1);

		assert(level_id[u] == level);
		assert(SR_rid[u] == R_end);
		++R_end;
		level_id[u] = n;
		degree[u] = 0;
		for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
		{
			ui w = edges[i];
			assert(SR[SR_rid[w]] == w);
			if (SR_rid[w] < R_end)
			{
				++degree[w];
				++degree[u];
			}
		}

		assert(level_id[u] > level);
		if (SR_rid[u] != S_end)
		{
			swap_pos(S_end, SR_rid[u]);
		}
		bool can_add = vc->verify_Sbundle_by_maxFlowAlg2(SR, S_end + 1, pstart, pend, edges, s);
		if (can_add)
		{
			assert(SR_rid[u] == S_end);
			++S_end;
			sbundle_binary_BB(S_end, R_end, level + 1);
			--S_end;
		}
	}

	bool greedily_add_vertices_to_S(ui &S_end, ui &R_end, ui level)
	{
		while (true)
		{
			ui *candidates = S2;
			ui candidates_n = 0;
			for (ui i = S_end; i < R_end; i++)
			{
				ui u = SR[i];
				if (R_end - degree[u] > s)
					continue;

				ui neighbors_n = 0, non_neighbors_n = 0;
				get_neighbors_and_non_neighbors(u, 0, R_end, level, neighbors_n, non_neighbors_n);
				assert(non_neighbors_n < s);
				bool OK = true;
				for (ui j = 0; j < non_neighbors_n; j++)
					if (R_end - degree[nonneighbors[j]] > s)
					{
						OK = false;
						break;
					}
				if (OK)
					candidates[candidates_n++] = u;
			}

			if (!candidates_n)
				break;

			while (candidates_n)
			{
				ui u = candidates[--candidates_n];
				assert(SR_rid[u] >= S_end);
				if (SR_rid[u] >= R_end)
					return false;

				bool can_add_u = false;
				if (S_end + 1 <= s)
				{
					can_add_u = true;
				}
				else
				{
					swap_pos(S_end, SR_rid[u]);
					can_add_u = vc->verify_Sbundle_by_maxFlowAlg2(SR, S_end + 1, pstart, pend, edges, s);
				}
				if (!can_add_u || !move_u_to_S_with_prune(u, S_end, R_end, level))
					return false;
			}
		}
		return true;
	}

	bool greedily_add_nonneighbors(ui *candidates, ui candidates_n, ui &S_end, ui &R_end, ui level)
	{
		while (candidates_n)
		{
			ui u = candidates[--candidates_n];
			assert(SR_rid[u] >= S_end);

			bool can_add_u = false;
			if (S_end + 1 <= s)
			{
				can_add_u = true;
			}
			else
			{
				swap_pos(S_end, SR_rid[u]);
				can_add_u = vc->verify_Sbundle_by_maxFlowAlg2(SR, S_end + 1, pstart, pend, edges, s);
			}
			if (SR_rid[u] >= R_end || !can_add_u || !move_u_to_S_with_prune(u, S_end, R_end, level))
				return false;
		}
		return true;
	}

	void get_neighbors_and_non_neighbors(ui u, ui idx_start, ui idx_end, ui level, ui &neighbors_n, ui &non_neighbors_n)
	{
		neighbors_n = non_neighbors_n = 0;
		for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
			if (SR_rid[edges[i]] < idx_end)
				vis[edges[i]] = 1;
		for (ui i = idx_start; i < idx_end; i++)
			if (SR[i] != u)
			{
				if (vis[SR[i]])
					neighbors[neighbors_n++] = SR[i];
				else
					nonneighbors[non_neighbors_n++] = SR[i];
			}
		for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
			vis[edges[i]] = 0;
	}

	bool move_u_to_S_with_prune(ui u, ui &S_end, ui &R_end, ui level)
	{
		assert(SR_rid[u] >= S_end && SR_rid[u] < R_end && SR[SR_rid[u]] == u);
		assert(degree_in_S[u] + s > S_end);
#ifndef NDEBUG
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SR[i]] + s >= S_end);
		for (ui i = S_end; i < R_end; i++)
			assert(degree_in_S[SR[i]] + s > S_end);
#endif
		if (SR_rid[u] != S_end)
			swap_pos(S_end, SR_rid[u]);
		++S_end;

		ui neighbors_n = 0, nonneighbors_n = 0;
		get_neighbors_and_non_neighbors(u, 0, R_end, level, neighbors_n, nonneighbors_n);
		assert(neighbors_n + nonneighbors_n == R_end - 1);
		for (ui i = 0; i < neighbors_n; i++)
			++degree_in_S[neighbors[i]];

		while (!Qv.empty())
			Qv.pop();

		for (ui i = 0; i < nonneighbors_n; i++)
		{
			ui v = nonneighbors[i];
			if (SR_rid[v] >= S_end)
			{
				if (level_id[v] == level)
					continue;
				if (S_end - degree_in_S[v] >= s || S_end - degree_in_S[u] == s)
				{
					level_id[v] = level;
					Qv.push(v);
				}
			}
			else if (S_end - degree_in_S[v] == s)
			{
				for (ept j = pstart[v]; j < pend[v] && level_id[edges[j]] >= level; j++)
					vis[edges[j]] = 1;
				for (ui j = S_end; j < R_end; j++)
					if (level_id[SR[j]] > level && !vis[SR[j]])
					{
						level_id[SR[j]] = level;
						Qv.push(SR[j]);
					}
				for (ept j = pstart[v]; j < pend[v] && level_id[edges[j]] >= level; j++)
					vis[edges[j]] = 0;
			}
		}
		return remove_vertices_with_prune(S_end, R_end, level);
	}

	bool remove_vertices_with_prune(ui S_end, ui &R_end, ui level)
	{
		while (true)
		{
			while (!Qv.empty())
			{
				ui u = Qv.front();
				Qv.pop();
				assert(SR[SR_rid[u]] == u);
				assert(SR_rid[u] >= S_end && SR_rid[u] < R_end);
				--R_end;
				swap_pos(SR_rid[u], R_end);

				bool terminate = false;
				ui neighbors_n = 0;

				for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
					if (SR_rid[edges[i]] < R_end)
					{

						ui w = edges[i];
						neighbors[neighbors_n++] = w;
						--degree[w];
						if (degree[w] + s <= best_solution_size)
						{
							if (SR_rid[w] < S_end)
								terminate = true;
							else if (level_id[w] > level)
							{
								level_id[w] = level;
								Qv.push(w);
							}
						}
					}
				if (terminate)
				{
					for (ui i = 0; i < neighbors_n; i++)
						++degree[neighbors[i]];
					level_id[u] = n;
					++R_end;
					return false;
				}
			}

#ifndef NDEBUG

			for (ui i = 0; i < R_end; i++)
			{
				assert(level_id[SR[i]] > level);
			}
#endif

			assert(Qv.empty());

			if (S_end >= 1)
			{
				std::vector<ui> deleted_vertices;
				delete_vertices_based_on_dist(S_end, R_end, s, deleted_vertices, level);
				if (!deleted_vertices.empty())
				{
					for (ui i = 0; i < deleted_vertices.size(); i++)
					{
						ui t_v = deleted_vertices[i];
						level_id[t_v] = level;
						assert(SR_rid[t_v] >= S_end && SR_rid[t_v] < R_end);
						Qv.push(t_v);
					}
				}
			}
			if (Qv.empty())
			{
				break;
			}
		}
		return true;
	}

	void restore_SR(ui &S_end, ui &R_end, ui old_S_end, ui old_R_end, ui level)
	{
		while (!Qv.empty())
		{
			ui u = Qv.front();
			Qv.pop();
			assert(level_id[u] == level);
			assert(SR_rid[u] < R_end);
			level_id[u] = n;
		}

		for (; R_end < old_R_end; R_end++)
		{
			ui u = SR[R_end];
			assert(level_id[u] == level && SR_rid[u] == R_end);
			level_id[u] = n;

			degree[u] = degree_in_S[u] = 0;
			for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
			{
				ui w = edges[i];
				assert(SR[SR_rid[w]] == w);
				if (SR_rid[w] < R_end)
				{
					++degree[w];
					++degree[u];
				}
				if (SR_rid[w] < S_end)
					++degree_in_S[u];
			}
		}

		for (; S_end > old_S_end; S_end--)
		{
			ui u = SR[S_end - 1];
			assert(SR_rid[u] == S_end - 1);

			for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
			{
				ui w = edges[i];
				if (SR_rid[w] < R_end)
					--degree_in_S[w];
			}
		}
	}

	void delete_vertices_based_on_dist(ui S_end, ui R_end, ui s, std::vector<ui> &results, ui level)
	{

		std::vector<ui> dist(R_end, n);
		std::queue<ui> qq;
		assert(S_end > 0 && s >= 2);
		ui tt_u = SR[S_end - 1];
		assert(SR_rid[tt_u] < S_end);
		dist[S_end - 1] = 0;
		qq.push(tt_u);
		while (!qq.empty())
		{
			ui current = qq.front();
			qq.pop();
			for (ui i = pstart[current]; i < pend[current] && level_id[edges[i]] >= level; ++i)
			{
				ui neighbor = edges[i];
				ui neighbor_id = SR_rid[neighbor];
				if (neighbor_id < R_end && dist[neighbor_id] == n)
				{
					dist[neighbor_id] = dist[SR_rid[current]] + 1;
					qq.push(neighbor);
				}
			}
		}

		for (ui i = S_end; i < R_end; i++)
		{
			assert(s >= 2);

			if (dist[i] > 2 + (s - 2) / (best_solution_size + 1 - s))
			{
				results.push_back(SR[i]);
			}
		}
	}
	bool remove_u_from_S_with_prune(ui &S_end, ui &R_end, ui level)
	{
		assert(S_end);
		ui u = SR[S_end - 1];
		--S_end;
		--R_end;
		swap_pos(S_end, R_end);
		level_id[u] = level;

		bool terminate = false;
		for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
			if (SR_rid[edges[i]] < R_end)
			{
				ui v = edges[i];
				--degree_in_S[v];
				--degree[v];
				if (degree[v] + s <= best_solution_size)
				{
					if (SR_rid[v] < S_end)
						terminate = true;
					else
					{
						assert(level_id[v] > level);
						level_id[v] = level;
						Qv.push(v);
					}
				}
			}
		if (terminate)
			return false;
		return true;
	}

	bool remove_u_from_C_with_prune(ui deleted_vertex, ui S_end, ui &R_end, ui level)
	{
		assert(S_end);
		ui u = deleted_vertex;
		--R_end;
		swap_pos(SR_rid[u], R_end);
		level_id[u] = level;

		bool terminate = false;
		for (ept i = pstart[u]; i < pend[u] && level_id[edges[i]] >= level; i++)
			if (SR_rid[edges[i]] < R_end)
			{
				ui v = edges[i];
				--degree[v];
				if (degree[v] + s <= best_solution_size)
				{
					if (SR_rid[v] < S_end)
						terminate = true;
					else
					{
						assert(level_id[v] > level);
						level_id[v] = level;
						Qv.push(v);
					}
				}
			}
		if (terminate)
			return false;
		return true;
	}

	bool collect_removable_vertices(ui S_end, ui R_end, ui level)
	{
		for (ui i = 0; i < S_end; i++)
			if (degree[SR[i]] + s <= best_solution_size)
				return false;

		for (ui i = S_end; i < R_end; i++)
			if (level_id[SR[i]] > level)
			{
				ui v = SR[i];
				if (S_end - degree_in_S[v] >= s || degree[v] + s <= best_solution_size)
				{
					level_id[v] = level;
					Qv.push(v);
					continue;
				}
			}

		return true;
	}

	void swap_pos(ui i, ui j)
	{
		std::swap(SR[i], SR[j]);
		SR_rid[SR[i]] = i;
		SR_rid[SR[j]] = j;
	}

	void swap_pos2(ui i, ui j)
	{
		std::swap(SR_remap[i], SR_remap[j]);
		SR_remap_rid[SR_remap[i]] = i;
		SR_remap_rid[SR_remap[j]] = j;
	}

	ui choose_branch_vertex(ui S_end, ui R_end, ui level)
	{
		ui *D = buf;
		ui D_n = 0;
		for (ui i = 0; i < R_end; i++)
			if (R_end - degree[SR[i]] > s)
				D[D_n++] = SR[i];
		assert(D_n != 0);

		ui min_degree_in_S = n;
		for (ui i = 0; i < D_n; i++)
			if (degree_in_S[D[i]] < min_degree_in_S)
				min_degree_in_S = degree_in_S[D[i]];
		ui u = n, min_degree = n;
		for (ui i = 0; i < D_n; i++)
			if (degree_in_S[D[i]] == min_degree_in_S && degree[D[i]] < min_degree)
			{
				min_degree = degree[D[i]];
				u = D[i];
			}
		assert(u != n);

		if (SR_rid[u] < S_end)
		{
			ui max_degree = 0, b = n;
			ui neighbors_n = 0, nonneighbors_n = 0;
			get_neighbors_and_non_neighbors(u, S_end, R_end, level, neighbors_n, nonneighbors_n);
			assert(nonneighbors_n);
			for (ui i = 0; i < nonneighbors_n; i++)
				if (degree[nonneighbors[i]] > max_degree)
				{
					max_degree = degree[nonneighbors[i]];
					b = nonneighbors[i];
				}
			return b;
		}
		else if (degree_in_S[u] < S_end || R_end - degree[u] > s + 1)
			return u;
		else
		{
			ui max_degree = 0, w = n;
			for (ui i = S_end; i < R_end; i++)
				if (degree[SR[i]] > max_degree)
				{
					max_degree = degree[SR[i]];
					w = SR[i];
				}
			if (degree[w] + 1 >= R_end)
			{
				printf("!!! WA degree[w]: %u, R_end: %u\n", degree[w], R_end);
			}
			assert(degree[w] + 1 < R_end);
			if (R_end - degree[w] == 2)
				return w;
			ui neighbors_n = 0, nonneighbors_n = 0;
			get_neighbors_and_non_neighbors(w, S_end, R_end, level, neighbors_n, nonneighbors_n);
			assert(nonneighbors_n);
			for (ui i = 0; i < nonneighbors_n; i++)
				if (R_end - degree[nonneighbors[i]] == s + 1)
					return nonneighbors[i];
		}

		printf("!!! WA in choose_branch_vertex\n");
		return n;
	}

	ui choose_branch_vertex_based_on_non_neighbors(ui S_end, ui R_end, ui vertex_in_S)
	{
		assert(SR_rid[vertex_in_S] < S_end);
		ui u = n, min_degree_in_S = n;
		for (ui i = S_end; i < R_end; i++)
		{
			ui v = SR[i];
			bool is_adj = false;
			for (ui j = pstart[vertex_in_S]; j < pend[vertex_in_S]; j++)
			{
				if (edges[j] == v)
				{
					is_adj = true;
					break;
				}
			}

			if (!is_adj)
			{
				if (degree_in_S[v] < min_degree_in_S)
				{
					u = v;
					min_degree_in_S = degree_in_S[v];
				}
				else if (degree_in_S[v] == min_degree_in_S)
				{
					if (degree[v] > degree[u])
					{
						u = v;
						min_degree_in_S = degree_in_S[v];
					}
				}
			}
		}

#ifndef NDEBUG
		bool is_adj = false;
		for (ui j = pstart[vertex_in_S]; j < pend[vertex_in_S]; j++)
		{
			if (edges[j] == u)
			{
				is_adj = true;
				break;
			}
		}
		assert(u != n && !is_adj);
#endif
		return u;
	}

	bool binary_find(ui u, ui w, ept begin, ept end, ui *edges)
	{
		if (begin >= end)
			return 0;

		ui idx;
		while (begin + 1 < end)
		{
			idx = begin + (end - begin) / 2;
			if (edges[idx] > w)
				end = idx;
			else
				begin = idx;
		}

		if (edges[begin] == w)
		{
			return 1;
		}

		return 0;
	}
};
#endif
