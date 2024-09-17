#ifndef _SBUNDLE_BB_MATRIX_
#define _SBUNDLE_BB_MATRIX_

#include "Utility.h"
#include "Timer.h"
#include "sbundle_tool.h"

class SBUNDLE_BB_matrix
{
private:
	long long n;
	ui original_R_end;

	char *matrix;
	long long matrix_size;
	std::vector<ui> temp_deleted;
	MaximumFlow *vc;
	long long branch_count;

	ui *degree;
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
	std::queue<ui> Qv;
	ui *level_id;
	ui max_level;

	std::vector<std::pair<ui, ui>> vp;
	std::vector<ui> non_adj;

public:
	SBUNDLE_BB_matrix()
	{
		n = 0;
		matrix = nullptr;
		matrix_size = 0;

		degree = degree_in_S = nullptr;

		best_solution = nullptr;
		s = best_solution_size = _UB_ = 0;

		neighbors = nonneighbors = nullptr;
		S2 = nullptr;

		SR = SR_rid = nullptr;
		level_id = nullptr;
		max_level = 0;
	}

	~SBUNDLE_BB_matrix()
	{
		if (matrix != NULL)
		{
			delete[] matrix;
			matrix = NULL;
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
		if (level_id != NULL)
		{
			delete[] level_id;
			level_id = NULL;
		}
	}

	void allocateMemory(ui n)
	{
		if (n <= 0)
			return;

		matrix_size = 1;
		matrix = new char[matrix_size];
		degree = new ui[n];
		degree_in_S = new ui[n];
		best_solution = new ui[n];
		SR = new ui[n];
		SR_rid = new ui[n];
		neighbors = new ui[n];
		nonneighbors = new ui[n];
		S2 = new ui[n];
		level_id = new ui[n];
	}

	void load_graph(ui _n, const std::vector<std::pair<ui, ui>> &vp)
	{
		n = _n;
		if (((long long)n) * n > matrix_size)
		{
			do
			{
				matrix_size *= 2;
			} while (((long long)n) * n > matrix_size);
			delete[] matrix;
			matrix = new char[matrix_size];
		}
		memset(matrix, 0, sizeof(char) * ((long long)n) * n);
		for (ui i = 0; i < n; i++)
			degree[i] = 0;
		for (ui i = 0; i < vp.size(); i++)
		{
			assert(vp[i].first >= 0 && vp[i].first < n && vp[i].second >= 0 && vp[i].second < n);
			ui a = vp[i].first, b = vp[i].second;
			degree[a]++;
			degree[b]++;
			if (matrix[a * n + b])
			{
				printf("Duplicate edge in load_graph()\n");
			}

			matrix[a * n + b] = matrix[b * n + a] = 1;
		}
	}

	void sBundle(ui S_, ui UB_, std::vector<ui> &sbundle, bool must_include_0, const int m, long long &branch_node)
	{
		s = S_;
		branch_count = 0;
		_UB_ = UB_;

		best_solution_size = sbundle.size();
		ui R_end;
		initialization(R_end, must_include_0, true);
		vc = new MaximumFlow();
		vc->prepare(n, m);

		if (R_end && best_solution_size < _UB_)
			splex_MultiBB(0, R_end, 1, must_include_0);

		branch_node += branch_count;

		if (vc != nullptr)
		{
			delete vc;
			vc = nullptr;
		}

		if (best_solution_size > sbundle.size())
		{
			sbundle.clear();
			for (int i = 0; i < best_solution_size; i++)
				sbundle.push_back(best_solution[i]);
		}
	}

	void initialization(ui &R_end, bool must_include_0, bool is_precise)
	{

		ui *peel_sequence = neighbors;
		ui *core = nonneighbors;
		ui *vis = SR;
		memset(vis, 0, sizeof(ui) * n);
		ui max_core = 0, idx = n;
		for (ui i = 0; i < n; i++)
		{
			ui u, min_degree = n;
			ui min2_degree = n;
			for (ui j = 0; j < n; j++)
			{
				if (!vis[j] && degree[j] < min_degree)
				{
					u = j;
					min2_degree = min_degree;
					min_degree = degree[j];
				}
				else if (!vis[j] && degree[j] < min2_degree)
				{
					min2_degree = degree[j];
				}
			}
			if (min_degree > max_core)
				max_core = min_degree;
			core[u] = max_core;
			peel_sequence[i] = u;
			vis[u] = 1;

			for (ui j = 0; j < n; j++)
				if (!vis[j] && matrix[u * n + j])
					--degree[j];
		}

		if (is_precise)
		{
			R_end = 0;
			for (ui i = 0; i < n; i++)
				SR_rid[i] = n;
			for (ui i = 0; i < n; i++)
				if (core[i] + s > best_solution_size)
				{
					SR[R_end] = i;
					SR_rid[i] = R_end;
					++R_end;
				}

			if ((must_include_0 && SR_rid[0] == n) || best_solution_size >= _UB_)
			{
				R_end = 0;
				return;
			}

			for (ui i = 0; i < R_end; i++)
			{
				ui u = SR[i];
				degree[u] = degree_in_S[u] = 0;
				for (ui j = 0; j < R_end; j++)
					if (matrix[u * n + SR[j]])
						++degree[u];
			}

			memset(level_id, 0, sizeof(ui) * n);
			for (ui i = 0; i < R_end; i++)
				level_id[SR[i]] = n;

			if (!Qv.empty())
				printf("!!! Something wrong. Qv must be empty in initialization\n");
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

	void splex_MultiBB(ui S_end, ui R_end, ui level, bool choose_zero)
	{
		if (S_end > best_solution_size)
			store_solution(S_end);
		if (R_end > best_solution_size && is_kplex(R_end) && vc->verify_SBundle_by_MaxFlowAlg0(SR, R_end, matrix, n, s, true))
			store_solution(R_end);
		if (R_end <= best_solution_size + 1 || best_solution_size >= _UB_)
			return;

		branch_count++;
#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < R_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d2;
			assert(d1 == degree_in_S[SR[i]]);
			assert(d2 == degree[SR[i]]);
		}
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SR[i]] + s >= S_end);
		for (ui i = S_end; i < R_end; i++)
			assert(degree_in_S[SR[i]] + s > S_end);
		for (ui i = 0; i < S_end; i++)
			if (degree_in_S[SR[i]] + s == S_end)
			{
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = S_end; j < R_end; j++)
					assert(t_matrix[SR[j]]);
			}
		for (ui i = 0; i < R_end; i++)
			assert(level_id[SR[i]] > level);
#endif

#ifndef NDEBUG

#endif
		ui old_S_end = S_end, old_R_end = R_end;
		assert(Qv.empty());

		if (choose_zero && SR_rid[0] < R_end && !move_u_to_S_with_prune(0, S_end, R_end, level))
		{
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
			return;
		}

		choose_zero = false;

		ui S2_n = 0;
		for (ui i = 0; i < S_end; i++)
			if (R_end - degree[SR[i]] > s)
				S2[S2_n++] = SR[i];

		if (R_end > best_solution_size && is_kplex(R_end) && vc->verify_SBundle_by_MaxFlowAlg0(SR, R_end, matrix, n, s, true))
			store_solution(R_end);
		if (R_end <= best_solution_size + 1 || best_solution_size >= _UB_)
		{
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
			return;
		}

#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < R_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d2;
			assert(d1 == degree_in_S[SR[i]]);
			assert(d2 == degree[SR[i]]);
		}
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SR[i]] + s >= S_end);
		for (ui i = S_end; i < R_end; i++)
			assert(degree_in_S[SR[i]] + s > S_end);
		for (ui i = 0; i < S_end; i++)
			if (degree_in_S[SR[i]] + s == S_end)
			{
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = S_end; j < R_end; j++)
					assert(t_matrix[SR[j]]);
			}
		for (ui i = 0; i < R_end; i++)
			assert(level_id[SR[i]] > level);
#endif
		ui u = n;
		ui min_deg = n;
		assert(!choose_zero);

		for (ui i = 0; i < R_end; i++)
			if (degree[SR[i]] < min_deg)
			{
				u = SR[i];
				min_deg = degree[SR[i]];
			}

		assert(u != n && SR_rid[u] < R_end && min_deg != n);

		if (min_deg >= R_end - s)
		{

			original_R_end = R_end;
#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < S_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d1;
				d2 = d1;
				for (ui j = S_end; j < R_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d2;
				assert(d1 == degree_in_S[SR[i]]);
				assert(d2 == degree[SR[i]]);
			}
#endif
			ui pre_S_end = S_end, pre_R_end = R_end;
			sbundle_binary_BB(S_end, R_end, level + 1);
			assert(S_end == pre_S_end && R_end == pre_R_end);
#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < S_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d1;
				d2 = d1;
				for (ui j = S_end; j < R_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d2;
				assert(d1 == degree_in_S[SR[i]]);
				assert(d2 == degree[SR[i]]);
			}
#endif
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
			return;
		}

		if (SR_rid[u] < S_end)
			u = choose_branch_vertex_based_on_non_neighbors(S_end, R_end, u);

		assert(SR_rid[u] >= S_end && SR_rid[u] < R_end);
		assert(degree[u] + s > best_solution_size && degree[u] + s > S_end);

		swap_pos(S_end, SR_rid[u]);
		bool can_add_u = vc->verify_SBundle_by_MaxFlowAlg0(SR, S_end + 1, matrix, n, s, false);

		if (can_add_u)
		{

			ui pre_best_solution_size = best_solution_size, t_old_S_end = S_end, t_old_R_end = R_end;
			assert(Qv.empty());

#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < S_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d1;
				d2 = d1;
				for (ui j = S_end; j < R_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d2;
				assert(d1 == degree_in_S[SR[i]]);
				assert(d2 == degree[SR[i]]);
			}
			for (ui i = 0; i < S_end; i++)
				assert(degree_in_S[SR[i]] + s >= S_end);
			for (ui i = S_end; i < R_end; i++)
				assert(degree_in_S[SR[i]] + s > S_end);
			for (ui i = 0; i < S_end; i++)
				if (degree_in_S[SR[i]] + s == S_end)
				{
					char *t_matrix = matrix + SR[i] * n;
					for (ui j = S_end; j < R_end; j++)
						assert(t_matrix[SR[j]]);
				}
			for (ui i = 0; i < R_end; i++)
				assert(level_id[SR[i]] > level);
#endif
			if (move_u_to_S_with_prune(u, S_end, R_end, level))
			{
#ifndef NDEBUG
				for (ui i = 0; i < R_end; i++)
				{
					ui d1 = 0, d2 = 0;
					for (ui j = 0; j < S_end; j++)
						if (matrix[SR[i] * n + SR[j]])
							++d1;
					d2 = d1;
					for (ui j = S_end; j < R_end; j++)
						if (matrix[SR[i] * n + SR[j]])
							++d2;
					assert(d1 == degree_in_S[SR[i]]);
					assert(d2 == degree[SR[i]]);
				}
				for (ui i = 0; i < S_end; i++)
					assert(degree_in_S[SR[i]] + s >= S_end);
				for (ui i = S_end; i < R_end; i++)
					assert(degree_in_S[SR[i]] + s > S_end);
				for (ui i = 0; i < S_end; i++)
					if (degree_in_S[SR[i]] + s == S_end)
					{
						char *t_matrix = matrix + SR[i] * n;
						for (ui j = S_end; j < R_end; j++)
							assert(t_matrix[SR[j]]);
					}
				for (ui i = 0; i < R_end; i++)
					assert(level_id[SR[i]] > level);
#endif
				splex_MultiBB(S_end, R_end, level + 1, false);
			}

			if (best_solution_size >= _UB_)
				return;
			assert(S_end == t_old_S_end + 1 && SR[S_end - 1] == u);
			restore_SR(S_end, R_end, S_end, t_old_R_end, level);

#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < S_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d1;
				d2 = d1;
				for (ui j = S_end; j < R_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d2;
				assert(d1 == degree_in_S[SR[i]]);
				assert(d2 == degree[SR[i]]);
			}
			for (ui i = 0; i < R_end; i++)
				assert(level_id[SR[i]] > level);
#endif

			assert(Qv.empty());
			bool succeed = remove_u_from_S_with_prune(S_end, R_end, level);
			if (succeed && best_solution_size > pre_best_solution_size)
				succeed = collect_removable_vertices(S_end, R_end, level);
			if (succeed)
				succeed = remove_vertices_and_prune(S_end, R_end, level);

			if (succeed)
				splex_MultiBB(S_end, R_end, level + 1, false);

			if (best_solution_size >= _UB_)
				return;
			assert(S_end >= old_S_end && R_end <= old_R_end);
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
#ifndef NDEBUG
			for (ui i = 0; i < R_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < S_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d1;
				d2 = d1;
				for (ui j = S_end; j < R_end; j++)
					if (matrix[SR[i] * n + SR[j]])
						++d2;
				assert(d1 == degree_in_S[SR[i]]);
				assert(d2 == degree[SR[i]]);
			}
			for (ui i = 0; i < S_end; i++)
				assert(degree_in_S[SR[i]] + s >= S_end);
			for (ui i = S_end; i < R_end; i++)
				assert(degree_in_S[SR[i]] + s > S_end);
			for (ui i = 0; i < S_end; i++)
				if (degree_in_S[SR[i]] + s == S_end)
				{
					char *t_matrix = matrix + SR[i] * n;
					for (ui j = S_end; j < R_end; j++)
						assert(t_matrix[SR[j]]);
				}
			for (ui i = 0; i < R_end; i++)
				assert(level_id[SR[i]] > level);
#endif
		}
		else
		{
			ui pre_best_solution_size = best_solution_size, t_old_S_end = S_end, t_old_R_end = R_end;
			assert(Qv.empty());
			bool succeed = remove_u_from_C_with_prune(u, S_end, R_end, level);
			if (succeed && best_solution_size > pre_best_solution_size)
				succeed = collect_removable_vertices(S_end, R_end, level);
			if (succeed)
				succeed = remove_vertices_and_prune(S_end, R_end, level);

			if (succeed)
				splex_MultiBB(S_end, R_end, level + 1, false);

			if (best_solution_size >= _UB_)
				return;
			assert(S_end >= old_S_end && R_end <= old_R_end);
			restore_SR(S_end, R_end, old_S_end, old_R_end, level);
		}
	}

	void sbundle_binary_BB(ui S_end, ui R_end, ui level)
	{

		if (S_end > best_solution_size)
			store_solution(S_end);
		if (R_end > best_solution_size && is_kplex(R_end) && vc->verify_SBundle_by_MaxFlowAlg0(SR, R_end, matrix, n, s, true))
			store_solution(R_end);
		if (R_end <= best_solution_size + 1)
			return;

		branch_count++;

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
		delfrC(u, R_end, level);
		sbundle_binary_BB(S_end, R_end, level + 1);
		addtoC(u, R_end, level);

		assert(level_id[u] > level);
		if (SR_rid[u] != S_end)
		{
			swap_pos(S_end, SR_rid[u]);
		}
		bool can_add = vc->verify_SBundle_by_MaxFlowAlg0(SR, S_end + 1, matrix, n, s, false);
		if (can_add)
		{
			addtoS(u, S_end, R_end);
			sbundle_binary_BB(S_end, R_end, level + 1);
			delfrS(u, S_end, R_end);
		}
	}

	bool canadd(ui u, ui S_end, bool precise)
	{
		char *t_matrix = matrix + u * n;

		if (precise)
		{
			if (SR_rid[u] != S_end)
			{
				swap_pos(S_end, SR_rid[u]);
			}
			return vc->verify_SBundle_by_MaxFlowAlg0(SR, S_end + 1, matrix, n, s, false);
		}

		return true;
	}

	void delfrC(ui u, ui &R_end, ui level)
	{
		assert(level_id[u] == n);
		level_id[u] = level;
		--R_end;

		swap_pos(R_end, SR_rid[u]);

		char *t_matrix = matrix + u * n;
		for (ui i = 0; i < R_end; i++)
			if (t_matrix[SR[i]])
			{
				assert(degree[SR[i]] > 0);
				--degree[SR[i]];
			}
	}

	void addtoC(ui u, ui &R_end, ui level)
	{

		assert(level_id[u] != n && level_id[u] == level);
		level_id[u] = n;
		assert(SR_rid[u] == R_end);
		SR_rid[u] = R_end;
		SR[R_end] = u;
		++R_end;

		char *t_matrix = matrix + u * n;
		degree[u] = 0;
		for (ui i = 0; i < R_end; i++)
		{
			if (t_matrix[SR[i]])
			{
				++degree[SR[i]];
				degree[u]++;
			}
		}
	}

	void addtoS(ui u, ui &S_end, ui R_end)
	{
		swap_pos(S_end, SR_rid[u]);
		++S_end;
	}

	void delfrS(ui u, ui &S_end, ui R_end)
	{
		--S_end;
		assert(S_end == SR_rid[u]);
	}

	void collect_removable_vertices_based_on_total_edges(ui S2_n, ui S_end, ui R_end, ui level)
	{
		vp.resize(R_end - S_end);
		ui max_nn = 0;
		for (ui i = S_end; i < R_end; i++)
		{
			ui nn = 0;
			if (S2_n != S_end)
			{
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = 0; j < S2_n; j++)
					if (!t_matrix[S2[j]])
						++nn;
			}
			else
				nn = S_end - degree_in_S[SR[i]];
			if (nn > max_nn)
				max_nn = nn;
			vp[i - S_end].first = SR[i];
			vp[i - S_end].second = nn;
		}
		ui *cnt = neighbors;
		for (ui i = 0; i <= max_nn; i++)
			cnt[i] = 0;
		for (ui i = 0; i < vp.size(); i++)
			++cnt[vp[i].second];
		for (ui i = 0; i < max_nn; i++)
			cnt[i + 1] += cnt[i];
		for (ui i = max_nn; i > 0; i--)
			cnt[i] = cnt[i - 1];
		cnt[0] = 0;
		ui *ids = nonneighbors;
		for (ui i = 0; i < vp.size(); i++)
			ids[cnt[vp[i].second]++] = i;

		ui total_support = 0;
		for (ui i = 0; i < S2_n; i++)
			total_support += s - S_end + degree_in_S[S2[i]];

		ui new_n = 0;
		while (!Qv.empty())
			Qv.pop();
		for (ui i = 0; i < vp.size(); i++)
		{
			ui idx = ids[i], v = vp[ids[i]].first;
			ui t_support = total_support - vp[idx].second;
			char *t_matrix = matrix + v * n;
			ui j = 0, v_support = s - 1 - S_end + degree_in_S[v], ub = S_end + 1;
			while (true)
			{
				if (j == new_n)
					j = i + 1;
				if (j >= vp.size() || ub > best_solution_size || ub + vp.size() - j <= best_solution_size)
					break;
				ui u = vp[ids[j]].first, nn = vp[ids[j]].second;
				if (t_support < nn)
					break;
				if (t_matrix[u])
				{
					t_support -= nn;
					++ub;
				}
				else if (v_support > 0)
				{
					--v_support;
					t_support -= nn;
					++ub;
				}
				++j;
			}
			if (ub <= best_solution_size)
			{
				level_id[v] = level;
				Qv.push(v);
			}
			else
				ids[new_n++] = ids[i];
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

				char *t_matrix = matrix + u * n;
				bool OK = true;
				for (ui j = 0; j < R_end; j++)
					if (j != i && !t_matrix[SR[j]] && R_end - degree[SR[j]] > s)
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

				if (!move_u_to_S_with_prune(u, S_end, R_end, level))
					return false;
			}
		}
		return true;
	}

	void check_degrees(ui S_end, ui R_end)
	{
		for (ui i = 0; i < R_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < R_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d2;
			assert(d1 == degree_in_S[SR[i]]);
			assert(d2 == degree[SR[i]]);
		}
	}

	bool greedily_add_nonneighbors(ui *candidates, ui candidates_n, ui &S_end, ui &R_end, ui level)
	{
		while (candidates_n)
		{
			ui u = candidates[--candidates_n];
			assert(SR_rid[u] >= S_end);
			if (SR_rid[u] >= R_end || !move_u_to_S_with_prune(u, S_end, R_end, level))
				return false;
		}
		return true;
	}

	void get_neighbors_and_nonneighbors(ui u, ui R_end, ui &neighbors_n, ui &nonneighbors_n)
	{
		neighbors_n = 0;
		nonneighbors_n = 0;
		char *t_matrix = matrix + u * n;
		for (ui i = 0; i < R_end; i++)
			if (SR[i] != u)
			{
				if (t_matrix[SR[i]])
					neighbors[neighbors_n++] = SR[i];
				else
					nonneighbors[nonneighbors_n++] = SR[i];
			}
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
		for (ui i = 0; i < R_end; i++)
			assert(level_id[SR[i]] > level);
#endif
		if (SR_rid[u] != S_end)
			swap_pos(S_end, SR_rid[u]);
		++S_end;

		ui neighbors_n = 0, nonneighbors_n = 0;
		get_neighbors_and_nonneighbors(u, R_end, neighbors_n, nonneighbors_n);
		assert(neighbors_n + nonneighbors_n == R_end - 1);
		for (ui i = 0; i < neighbors_n; i++)
			++degree_in_S[neighbors[i]];
#ifndef NDEBUG
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SR[i]] + s >= S_end);
#endif

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
				char *tt_matrix = matrix + v * n;
				for (ui j = S_end; j < R_end; j++)
					if (level_id[SR[j]] > level && !tt_matrix[SR[j]])
					{
						level_id[SR[j]] = level;
						Qv.push(SR[j]);
					}
			}
		}

#ifndef NDEBUG
		for (ui i = 0; i < S_end; i++)
			if (degree_in_S[SR[i]] + s == S_end)
			{
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = S_end; j < R_end; j++)
					assert(level_id[SR[j]] == level || t_matrix[SR[j]]);
			}
#endif

		return remove_vertices_and_prune(S_end, R_end, level);
	}

	bool remove_vertices_and_prune(ui S_end, ui &R_end, ui level)
	{
		while (true)
		{
			while (!Qv.empty())
			{
				ui u = Qv.front();
				Qv.pop();
				assert(SR[SR_rid[u]] == u && SR_rid[u] >= S_end && SR_rid[u] < R_end);
				--R_end;
				swap_pos(SR_rid[u], R_end);

				bool terminate = false;
				ui neighbors_n = 0;
				char *t_matrix = matrix + u * n;
				for (ui i = 0; i < R_end; i++)
					if (t_matrix[SR[i]])
					{
						ui w = SR[i];
						neighbors[neighbors_n++] = w;
						--degree[w];
						if (degree[w] + s <= best_solution_size)
						{
							if (i < S_end)
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

			if (Qv.empty())
			{
				break;
			}
		}

		assert(Qv.empty());
		return true;
	}

	void restore_SR(ui &S_end, ui &R_end, ui old_S_end, ui old_R_end, ui level)
	{
#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < R_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d2;
			assert(d1 == degree_in_S[SR[i]]);
			assert(d2 == degree[SR[i]]);
		}
#endif

		while (!Qv.empty())
		{
			ui u = Qv.front();
			Qv.pop();
			assert(level_id[u] == level);
			assert(SR_rid[u] < R_end);
			level_id[u] = n;
		}

#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			if (level_id[SR[i]] <= level)
				printf("level_id[%u] = %u, level = %u, n = %u\n", SR[i], level_id[SR[i]], level, n);
			assert(level_id[SR[i]] > level);
		}
#endif

		for (; R_end < old_R_end; R_end++)
		{
			ui u = SR[R_end];
			assert(level_id[u] == level && SR_rid[u] == R_end);
			level_id[u] = n;
			ui neighbors_n = 0;
			char *t_matrix = matrix + u * n;
			degree[u] = degree_in_S[u] = 0;
			for (ui i = 0; i < R_end; i++)
				if (t_matrix[SR[i]])
				{
					ui w = SR[i];
					neighbors[neighbors_n++] = w;
					++degree[w];
					++degree[u];
					if (i < S_end)
						++degree_in_S[u];
				}
		}

		for (; S_end > old_S_end; S_end--)
		{
			ui u = SR[S_end - 1];
			assert(SR_rid[u] == S_end - 1);

			ui neighbors_n = 0;
			char *t_matrix = matrix + u * n;
			for (ui i = 0; i < R_end; i++)
				if (t_matrix[SR[i]])
				{
					ui w = SR[i];
					neighbors[neighbors_n++] = w;
					--degree_in_S[w];
				}
		}

#ifndef NDEBUG
		for (ui i = 0; i < R_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < R_end; j++)
				if (matrix[SR[i] * n + SR[j]])
					++d2;
			assert(d1 == degree_in_S[SR[i]]);
			assert(d2 == degree[SR[i]]);
		}
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SR[i]] + s >= S_end);
		for (ui i = S_end; i < R_end; i++)
			assert(old_S_end == S_end || degree_in_S[SR[i]] + s > S_end);
		for (ui i = 0; i < S_end; i++)
			if (degree_in_S[SR[i]] + s == S_end)
			{
				char *t_matrix = matrix + SR[i] * n;
				for (ui j = S_end; j < R_end; j++)
					assert(old_S_end == S_end || t_matrix[SR[j]]);
			}
		for (ui i = 0; i < R_end; i++)
			assert(level_id[SR[i]] > level);
#endif
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
		ui neighbors_n = 0;
		char *t_matrix = matrix + u * n;
		for (ui i = 0; i < R_end; i++)
			if (t_matrix[SR[i]])
				neighbors[neighbors_n++] = SR[i];
		for (ui i = 0; i < neighbors_n; i++)
		{
			ui v = neighbors[i];
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
		ui neighbors_n = 0;
		char *t_matrix = matrix + u * n;
		for (ui i = 0; i < R_end; i++)
			if (t_matrix[SR[i]])
				neighbors[neighbors_n++] = SR[i];
		for (ui i = 0; i < neighbors_n; i++)
		{
			ui v = neighbors[i];
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
				char *t_matrix = matrix + v * n;
				for (ui j = 0; j < S_end; j++)
				{
					if (S_end - degree_in_S[SR[j]] == s && !t_matrix[SR[j]])

					{
						level_id[v] = level;
						Qv.push(v);
						break;
					}
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

	ui choose_branch_vertex(ui S_end, ui R_end)
	{
		ui *D = neighbors;
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
			char *t_matrix = matrix + u * n;
			for (ui i = S_end; i < R_end; i++)
				if (!t_matrix[SR[i]] && degree[SR[i]] > max_degree)
				{
					max_degree = degree[SR[i]];
					b = SR[i];
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
			char *t_matrix = matrix + w * n;
			for (ui i = S_end; i < R_end; i++)
				if (!t_matrix[SR[i]] && R_end - degree[SR[i]] == s + 1)
					return SR[i];
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
			if (!matrix[vertex_in_S * n + v])
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

		assert(u != n && !matrix[vertex_in_S * n + u]);
		return u;
	}
};
#endif
