/*
 * CTPrune.h
 *
 *  Created on: 31 May 2023
 *      Author: ljchang@outlook.com
 */

#ifndef CTPRUNE_H_
#define CTPRUNE_H_

#include <cstring>
#include <algorithm>

using ui = unsigned int;  // vertex type
using ept = unsigned int; // edge pointer type; unsigned int can be used to process upto two billion undirected edges

namespace CTPrune
{
	// orient graph. 注意，这里的做法就相当于是将之前无向图的m条边转化成了有向图的m/2条边。
	void orient_graph(ui n, ept m, ui *peel_sequence, ept *pstart, ept *pend, ui *edges, ui *rid)
	{
		// 注意，经过之前的调整，peel_sequence数组已经装好了正确且有效的顶点集合，并且索引是从0到cnt
		for (ui i = 0; i < n; i++)
			rid[peel_sequence[i]] = i; // 这里改变了rid映射数组的值，使得rid数组和peel_sequence中的元素一一对应
		for (ui i = 0; i < n; i++)
		{
			ept &end = pend[i] = pstart[i]; // 这里的&end是pend[i]的一个别名，主要是为了写起来方便
			for (ept j = pstart[i]; j < pstart[i + 1]; j++)
				if (rid[edges[j]] > rid[i])
					edges[end++] = edges[j]; // 我有点好奇，如果这里只保留有向边的话，即（3,4), (3,5)这种边，那假如之前的顶点3有一条边是(3,2)的话，按照这种有向的表示方式不是那条边直接就没了？不过我感觉还是可以通过pstart+1来判断
		}
	}

	// oriented triangle counting. 在常老师这个版本的实现中，三角形的计数变得格外的容易
	void oriented_triangle_counting(ui n, ept m, ept *pstart, ept *pend, ui *edges, ui *tri_cnt, ui *adj)
	{
		memset(adj, 0, sizeof(ui) * n); // adj数组用来暂时标记某个顶点u的所有邻居从而可以实现在O(1)的时间复杂度内判断另外一个顶点v是不是顶点u的邻居。这种设计在寻找共同邻居的时候特别有用
		// long long cnt = 0;
		memset(tri_cnt, 0, sizeof(ui) * m);
		for (ui u = 0; u < n; u++)
		{
			for (ept j = pstart[u]; j < pend[u]; j++)
				adj[edges[j]] = j + 1; // 先把顶点u的所有邻居都进行一个标记（仅在当前轮次的循环中）

			for (ept j = pstart[u]; j < pend[u]; j++)
			{
				ui v = edges[j];						  // 然后找到u的所有1跳邻居顶点v
				for (ept k = pstart[v]; k < pend[v]; k++) // 接着寻找v的所有1跳邻居顶点k。20240521：注意，之所以可以这么统计，是因为前面使用了orient_graph对边进行了有向化了。
					if (adj[edges[k]])					  // 如果顶点k也是之前我们所标记的顶点u的邻居顶点，那么就说明这时候找到了一个由u, v, k组成的三角形了
					{									  // 因此更新下面三者的计数
						++tri_cnt[j];
						++tri_cnt[k];
						++tri_cnt[adj[edges[k]] - 1];
						// ++ cnt;
					}
			}

			for (ept j = pstart[u]; j < pend[u]; j++)
				adj[edges[j]] = 0; // 清除在当前轮次中顶点u的所有邻居的标记
		}
	}

	// 该函数改变了图中的顶点数n和边数m
	bool remove_and_shrink_oriented_tri(ui &n, ept &m, ui triangle_threshold, ui *out_mapping, ui *peel_sequence, ept *pstart, ept *pend, ui *edges, ui *tri_cnt, ui *rid, ui *degree)
	{
		// 在该函数内，又对out_mapping进行了重新整理，即只保留度数大于threshold的顶点集合。20240502：所以说out_mapping中的顶点编号始终对应着最原始的顶点编号。
		for (ui i = 0; i < n; i++)
			degree[i] = pstart[i + 1] - pstart[i];
		ept removed_edges = 0; // 用来统计所删除的边的数量
		for (ui i = 0; i < n; i++)
			for (ui j = pstart[i]; j < pend[i]; j++)
				if (tri_cnt[j] < triangle_threshold)
				{
					--degree[i];
					--degree[edges[j]];
					++removed_edges;
				}

		if (removed_edges <= m / 4)
			return false; // 如果在这一次的移除过程中减少的边的数量少于当前边的总数的1/4，那么就返回false，用来告诉后面不用再继续调用这个函数了。应该是出于性能考虑，防止外面的while循环执行时间过长

		ui cnt = 0;
		for (ui i = 0; i < n; i++)
			if (degree[i] > 0) // 重新统计当前图中度数符合要求的顶点
			{
				out_mapping[cnt] = out_mapping[i]; // 这里采用的做法是直接把之前计算好的值挪过来，从而使得映射一直有效。所以说，out_mapping里面一直装的就是当前图中有效顶点的最原始的编号
				rid[i] = cnt++;					   // 这里rid的角色开始发生变化，rid用来反映某个顶点在out_mapping数组中的编号
			}

		ui t_cnt = 0;
		for (ui i = 0; i < n; i++)
			if (degree[peel_sequence[i]] > 0)
				peel_sequence[t_cnt++] = rid[peel_sequence[i]]; // 同时也维护peel_sequence中的顶点集合。经过该操作，peel_sequence中存储的其实就是out_mapping中的degeneracy ordering

		// 下面用来删除边的方式非常的巧妙，直接通过调整edges数组中的指针的位置（仅统计没有被删除过的边的位置）就表示把边已经删掉了
		ui pos = 0; // 在常老师的代码中，pos是用来跟踪边的变量，而cnt则是用来跟踪顶点的变量
		cnt = 0;
		for (ui i = 0; i < n; i++)
			if (degree[i] > 0)
			{
				ui start = pstart[i];
				pstart[cnt] = pos;
				for (ui j = start; j < pend[i]; j++)
					if (tri_cnt[j] >= triangle_threshold)
						edges[pos++] = rid[edges[j]]; // 慢慢理清楚了这里面的关系了。rid中装的是顶点在out_mapping中的索引，而out_mapping中装的则是最原始的顶点的索引。环环相扣
				pend[cnt] = pos;
				pos = degree[i] + pstart[cnt]; // 作为下一个顶点的pstart的值
				++cnt;
			}
		pstart[cnt] = m = pos; // 此处的边数和顶点数量再次发生了变化
		n = cnt;

		return true;
	}

	// reorganize the adjacency lists and sort each adjacency list to be in increasing order
	// 这个函数的作用是重新组织邻接表。经过该函数之后，pstart[i]和pstart[i+1]之间就是存放的顶点i当前所有有效邻居顶点
	void reorganize_oriented_graph(ui n, ui *tri_cnt, ept *pstart, ept *pend, ept *pend2, ui *edges, ept *edges_pointer, ui *buf)
	{
		for (ui i = 0; i < n; i++)
			pend2[i] = pend[i];
		for (ui i = 0; i < n; i++)
		{
			for (ept j = pstart[i]; j < pend[i]; j++)
			{
				ept &k = pend2[edges[j]]; // 一条边的一个端点
				edges[k] = i;			  // 该条边的另外一个端点，也就是i本身
				tri_cnt[k] = tri_cnt[j];
				++k; // 注意这里的k是对于pend2的一个引用，所以说这里其实是改变的pend2的值
			}
		}
		// 之所以定义一个pend2的原因是上面的代码中涉及到了修改pend2，而这个修改不方便直接改变原始的pend数组。

		for (ui i = 0; i < n; i++)
		{
			pend2[i] = pend[i]; // 由于上面的代码中改变了pend2，这里重新把正确的pend赋值过去
			pend[i] = pstart[i];
		}
		for (ui i = 0; i < n; i++)
		{
			for (ept j = pend2[i]; j < pstart[i + 1]; j++)
			{
				ept &k = pend[edges[j]]; // 注意，这里修改的就是pend数组了
				edges[k] = i;
				tri_cnt[k] = tri_cnt[j];
				edges_pointer[k] = j; // 编号为k的边指向j
				edges_pointer[j] = k; // 编号为j的边指向k
				++k;
			}
		}

		// 下面的代码就是重新组织edges、tri_cnt和edges_pointer数组
		ept *pointer = pend2;
		for (ui i = 0; i < n; i++)
		{
			if (pend[i] == pstart[i] || pend[i] == pstart[i + 1]) // 这里在常老师之前的kpelxS的代码中出现过
				continue;
			ept j = pstart[i], k = pend[i], pos = 0;
			while (j < pend[i] || k < pstart[i + 1])
			{
				if (k >= pstart[i + 1] || (j < pend[i] && edges[j] < edges[k]))
				{
					buf[pos] = edges[j];
					pointer[pos++] = edges_pointer[j++];
				}
				else
				{
					buf[pos] = edges[k];
					pointer[pos++] = edges_pointer[k++];
				}
			}
			for (ept j = 0; j < pos; j++)
			{
				ept idx = pstart[i] + j, k = pointer[j];
				edges[idx] = buf[j];
				tri_cnt[idx] = tri_cnt[k];
				edges_pointer[idx] = k; // 所以从这里可以看出来edges_pointer就是存储的每一条边的另外一个顶点。更具体来说，一共有m条边，这m条边的起点就是edges_pointer[idx]中的下标idx，而终点就是k。
				edges_pointer[k] = idx;
			}
		}
	}

	bool find(ui w, ept b, ept e, char *deleted, ept &idx, ui *edges) // 我感觉这是常老师在代码里面写的基于邻接表形式的二分查找函数。我进一步的确定了，感觉就是的。参考常老师对这个函数的调用，可以推测出形式非常像
	{
		if (b >= e)
			return false;

		while (b + 1 < e)
		{
			idx = b + (e - b) / 2;
			if (edges[idx] > w)
				e = idx;
			else
				b = idx;
		}

		idx = b;
		if (edges[idx] == w && !deleted[idx])
			return true;
		return false;
	}

	// 该函数的作用就是把pend[u]设置成pstart[u+1]，同时把对应的一些变量都顺移过去
	void compact_neighbors(ui u, ui *tri_cnt, ui *edges_pointer, char *deleted, ept *pstart, ept *pend, ui *edges)
	{
		ui end = pstart[u];
		for (ui i = pstart[u]; i < pend[u]; i++)
			if (!deleted[i]) // 如果该条边没有被删除
			{
				edges[end] = edges[i];
				tri_cnt[end] = tri_cnt[i];
				edges_pointer[end] = edges_pointer[i];	 // 这里表示一条有向边（end, edges_pointer[i])
				edges_pointer[edges_pointer[end]] = end; // 这里表示上面反向的一条有向边(edges_pointer[end], end)
				deleted[end] = 0;						 // 表示edges数组中的end位置没有被删除
				++end;
			}
		pend[u] = end; // 一般来说pstart都是不会变的，而只改变pend的指向
	}

	// 这里和之前常老师的代码不太一样了，这里是单独的把truss的缩减打包成了一个函数
	void truss_peeling(ui degree_threshold, ui *Qv, ui &Qv_n, ui triangle_threshold, ui *Qe, ept Qe_n, ui *tri_cnt, ept *edges_pointer, char *deleted, ui *degree, ept *pstart, ept *pend, ui *edges)
	{
		for (ept j = 0; j < Qe_n; j += 2)
		{
			ui u = Qe[j], v = Qe[j + 1], idx;
			find(v, pstart[u], pend[u], deleted, idx, edges); // 20240521：没想到常老师在代码里面也用到了二分查找的思想

			ui tri_n = tri_cnt[idx];
			deleted[idx] = deleted[edges_pointer[idx]] = 1;
			if (Qv != nullptr)
			{
				if (degree[u] == degree_threshold)
					Qv[Qv_n++] = u;
				if (degree[v] == degree_threshold)
					Qv[Qv_n++] = v;
			}
			--degree[u];
			--degree[v];
			if (pend[u] - pstart[u] > degree[u] * 2)
				compact_neighbors(u, tri_cnt, edges_pointer, deleted, pstart, pend, edges);
			if (pend[v] - pstart[v] > degree[v] * 2)
				compact_neighbors(v, tri_cnt, edges_pointer, deleted, pstart, pend, edges);
			if (pend[u] - pstart[u] < pend[v] - pstart[v])
				std::swap(u, v);

			if (pend[u] - pstart[u] > (pend[v] - pstart[v]) * 2)
			{ // binary search
				for (ept k = pstart[v]; k < pend[v]; k++)
					if (!deleted[k])
					{
						if (tri_n && find(edges[k], pstart[u], pend[u], deleted, idx, edges))
						{
							--tri_n;
							--tri_cnt[edges_pointer[idx]];
							if ((tri_cnt[idx]--) == triangle_threshold)
							{
								Qe[Qe_n++] = u;
								Qe[Qe_n++] = edges[idx];
							}
							--tri_cnt[edges_pointer[k]];
							if ((tri_cnt[k]--) == triangle_threshold)
							{
								Qe[Qe_n++] = v;
								Qe[Qe_n++] = edges[k];
							}
						}
					}
			}
			else
			{ // sorted_merge
				ept ii = pstart[u], jj = pstart[v];
				while (ii < pend[u] && jj < pend[v])
				{
					if (edges[ii] == edges[jj])
					{
						if (!deleted[ii] && !deleted[jj])
						{
							--tri_n;
							--tri_cnt[edges_pointer[ii]];
							if ((tri_cnt[ii]--) == triangle_threshold)
							{
								Qe[Qe_n++] = u;
								Qe[Qe_n++] = edges[ii];
							}
							--tri_cnt[edges_pointer[jj]];
							if ((tri_cnt[jj]--) == triangle_threshold)
							{
								Qe[Qe_n++] = v;
								Qe[Qe_n++] = edges[jj];
							}
						}

						++ii;
						++jj;
					}
					else if (edges[ii] < edges[jj])
						++ii;
					else
						++jj;
				}
			}
		}
	}

	void core_truss_copeeling(ui n, ui degree_threshold, ui triangle_threshold, ui *Qv, ui *Qe, ui *tri_cnt, ept *edges_pointer, char *deleted, char *exists, ui *degree, ept *pstart, ept *pend, ui *edges)
	{
		ept Qe_n = 0;
		for (ui i = 0; i < n; i++)
			for (ept j = pstart[i]; j < pend[i]; j++)				 // 所以说其实j就是代表着一个唯一的边的编号
				if (tri_cnt[j] < triangle_threshold && edges[j] > i) // 先判断有没有能够被删除的边
				{
					Qe[Qe_n++] = i; // 每次往边删除队列里面放入顶点的时候也是把编号小的顶点放在前面，编号大的放在后面
					Qe[Qe_n++] = edges[j];
				}
		ui Qv_n = 0;
		for (ui i = 0; i < n; i++)
			if (degree[i] < degree_threshold) // 再判断有没有能够删除的顶点
				Qv[Qv_n++] = i;
		while (Qe_n || Qv_n) // 这里又是在kplexS中比较熟悉的CTCP的流程了
		{
			while (Qe_n == 0 && Qv_n != 0)
			{
				ui u = Qv[--Qv_n]; // delete u from the graph due to have a degree < degree_threshold
				if (degree[u] == 0)
					continue;

				if (pend[u] - pstart[u] != degree[u])
					compact_neighbors(u, tri_cnt, edges_pointer, deleted, pstart, pend, edges); // 如果顶点u的度数不等于pstart和pend的跨度，那么就调整pend的位置
				assert(pend[u] - pstart[u] == degree[u]);

				for (ept i = pstart[u]; i < pend[u]; i++)
					deleted[i] = deleted[edges_pointer[i]] = exists[edges[i]] = 1; // 这里exists设置为1就是为了在下面的if语句中作为一个当前顶点存在的标记
				degree[u] = 0;

				for (ept i = pstart[u]; i < pend[u]; i++)
				{
					ui v = edges[i];
					if ((degree[v]--) == degree_threshold)
						Qv[Qv_n++] = v;
					if (pend[v] - pstart[v] > degree[v] * 2)										// 这里是指如果顶点v的邻接表（即通过pstart和pend共同表示）已经空出了很多位置了（通过实际度数degree来辅助判断），那么就优化邻居的分布，从而遍历效率更高
						compact_neighbors(v, tri_cnt, edges_pointer, deleted, pstart, pend, edges); // 似乎在很多地方都看到了这种优化技巧，这么做的目的是为了能够在下面代码的对顶点v的邻居遍历更加高效

					for (ept j = pstart[v]; j < pend[v]; j++)
						if (!deleted[j] && edges[j] > v && exists[edges[j]])
						{
							--tri_cnt[edges_pointer[j]];
							if ((tri_cnt[j]--) == triangle_threshold)
							{
								Qe[Qe_n++] = v, Qe[Qe_n++] = edges[j];
							}
						}
				}

				for (ept i = pstart[u]; i < pend[u]; i++)
					exists[edges[i]] = 0;
			}
			truss_peeling(degree_threshold, Qv, Qv_n, triangle_threshold, Qe, Qe_n, tri_cnt, edges_pointer, deleted, degree, pstart, pend, edges);
			Qe_n = 0;
		}
	}

	// reduce the graph to its maximal subgraph with and minimum edge triangle count at least triangle_threshold
	void truss_pruning(ui &n, ept &m, ui triangle_threshold, ui *peel_sequence, ui *out_mapping, ui *rid, ept *pstart, ui *edges, ui *degree, bool output)
	{
		if (triangle_threshold == 0)
		{
			printf("!!! Triangle_threshold is 0\n");
			return;
		}

		ept *pend = new ui[n + 1];
		orient_graph(n, m, peel_sequence, pstart, pend, edges, rid);
		ui *tri_cnt = new ui[m];
		oriented_triangle_counting(n, m, pstart, pend, edges, tri_cnt, rid);

		while (n && remove_and_shrink_oriented_tri(n, m, triangle_threshold, out_mapping, peel_sequence, pstart, pend, edges, tri_cnt, rid, degree))
		{
			oriented_triangle_counting(n, m, pstart, pend, edges, tri_cnt, rid);
		}

		ept *pend_buf = new ept[n + 1];
		ept *edges_pointer = new ept[m];
		reorganize_oriented_graph(n, tri_cnt, pstart, pend, pend_buf, edges, edges_pointer, rid);
		delete[] pend_buf;
		pend_buf = nullptr;

		for (ui i = 0; i < n; i++)
		{
			pend[i] = pstart[i + 1];
			degree[i] = pstart[i + 1] - pstart[i];
		}
		ui *Qe = new ui[m];
		char *deleted = new char[m];
		memset(deleted, 0, sizeof(char) * m);
		ept Qe_n = 0;
		for (ui i = 0; i < n; i++)
			for (ui j = pstart[i]; j < pend[i]; j++)
				if (tri_cnt[j] < triangle_threshold && edges[j] > i)
				{
					Qe[Qe_n++] = i;
					Qe[Qe_n++] = edges[j];
				}
		truss_peeling(0, nullptr, Qe_n, triangle_threshold, Qe, Qe_n, tri_cnt, edges_pointer, deleted, degree, pstart, pend, edges);

		ui cnt = 0;
		for (ui i = 0; i < n; i++)
			if (degree[i] > 0)
			{
				out_mapping[cnt] = out_mapping[i];
				rid[i] = cnt++;
			}
		ui t_cnt = 0;
		for (ui i = 0; i < n; i++)
			if (degree[peel_sequence[i]] > 0)
				peel_sequence[t_cnt++] = rid[peel_sequence[i]];
		ui pos = 0;
		cnt = 0;
		for (ui i = 0; i < n; i++)
			if (degree[i] > 0)
			{
				ui start = pstart[i];
				pstart[cnt] = pos;
				for (ui j = start; j < pend[i]; j++)
					if (!deleted[j])
						edges[pos++] = rid[edges[j]];
				++cnt;
			}
		pstart[cnt] = m = pos;
		n = cnt;

		delete[] Qe;
		delete[] deleted;
		delete[] edges_pointer;
		delete[] pend;
		delete[] tri_cnt;

		if (output)
			printf("*** After truss_pruning: n = %u, m = %lu (undirected)\n", n, m / 2);
	}

	// reduce the graph to its maximal subgraph with minimum degree at least degree_threshold and minimum edge triangle count at least triangle_threshold
	void core_truss_copruning(ui &n, ept &m, ui degree_threshold, ui triangle_threshold, ui *peel_sequence, ui *out_mapping, ui *rid, ept *pstart, ui *edges, ui *degree, bool output)
	{
		if (triangle_threshold == 0)
		{
			printf("!!! Triangle_threshold is 0\n");
			return;
		}
		
		if (degree_threshold <= triangle_threshold + 1)
		{
			printf("!!! Degree_threshold <= triangle_threshold + 1, please invoke truss_pruning\n"); // 如果当前的度数下界仅仅比truss的下界多了1，那么就应该直接调用truss-pruning，而不是这个CTCP
			return;
		}

		ept *pend = new ui[n + 1];									 // 在这里首次出现pend，用来表示一个顶点的邻接表的结束位置
		orient_graph(n, m, peel_sequence, pstart, pend, edges, rid); // 只有在先把图的邻接表进行有向排列之后，下面的三角形计数才避免了重复计数的问题。这是一种巧妙的做法。
		ui *tri_cnt = new ui[m];									 // 用来存储每一条边所参与构成的三角形计数
		oriented_triangle_counting(n, m, pstart, pend, edges, tri_cnt, rid);

		while (n && remove_and_shrink_oriented_tri(n, m, triangle_threshold, out_mapping, peel_sequence, pstart, pend, edges, tri_cnt, rid, degree))
		{
			oriented_triangle_counting(n, m, pstart, pend, edges, tri_cnt, rid);
		}

		ept *pend_buf = new ept[n + 1];
		ept *edges_pointer = new ept[m];
		reorganize_oriented_graph(n, tri_cnt, pstart, pend, pend_buf, edges, edges_pointer, rid); // 这里主要是把edges_pointer里面的值填好
		delete[] pend_buf;
		pend_buf = nullptr;

		for (ui i = 0; i < n; i++)
		{
			pend[i] = pstart[i + 1]; // 因为前面经过reorganize_oriented_graph函数，整个邻接表变得非常的紧凑，所以这里可以把pend的值设置成下一个顶点的pstart的位置
			degree[i] = pstart[i + 1] - pstart[i];
		}
		ui *Qe = new ui[m];
		ui *Qv = new ui[n];
		// 20240425：注意！！！原来这个deleted数组仅仅只是一个局部数组！即在CTCP中使用和被释放的！
		char *deleted = new char[m]; // 果然，deleted数组是用来表示哪些边被删除与否
		memset(deleted, 0, sizeof(char) * m);
		char *exists = new char[n];
		memset(exists, 0, sizeof(char) * n);
		core_truss_copeeling(n, degree_threshold, triangle_threshold, Qv, Qe, tri_cnt, edges_pointer, deleted, exists, degree, pstart, pend, edges);

		ui cnt = 0;
		for (ui i = 0; i < n; i++)
			if (degree[i] > 0)
			{
				out_mapping[cnt] = out_mapping[i]; // 这里又对当前图中存在的顶点进行了删除，即只保留那些degree[i] > 0的顶点
				rid[i] = cnt++;					   // 每一次调整完out_mapping之后，也会相应的调整rid数组，从而使得rid数组完整的保留了out_mapping中顶点在原始图中的映射
												   // 20240502：注意，out_mapping数组中始终对应着最原始顶点的编号。
			}
		ui t_cnt = 0;
		for (ui i = 0; i < n; i++)
			if (degree[peel_sequence[i]] > 0)
				peel_sequence[t_cnt++] = rid[peel_sequence[i]]; // 前面整理完了out_mapping和rid数组，这里又对peel_sequence数组进行整理。peel_sequence中装的就是out_mapping中顶点的退化序。这个序是不会变的。
		ui pos = 0;
		cnt = 0;
		for (ui i = 0; i < n; i++) // 下面开始整理pstart数组，即如果某些边要是被删除了，那么就把它们在edges数组中所对应的顶点也挪走。不过需要注意，这里只维护了pstart！没有维护pend！
			if (degree[i] > 0)
			{
				ui start = pstart[i];
				pstart[cnt] = pos;
				for (ui j = start; j < pend[i]; j++)
					if (!deleted[j])
						edges[pos++] = rid[edges[j]];
				++cnt;
			}
		pstart[cnt] = m = pos; // 由于该函数删除了一些顶点和边，所以这里更新边和顶点的数量
		n = cnt;

		delete[] Qv;
		delete[] Qe;
		delete[] deleted;
		delete[] edges_pointer;
		delete[] pend;
		delete[] tri_cnt;

		if (output)
			printf("*** After core_truss_copruning: n = %u, m = %lu (undirected)\n", n, m / 2);
	}

	void check_core_pruning(ui n, ui m, ui degree_threshold, ept *pstart, ui *edges)
	{
		for (ui i = 0; i < n; i++)
		{
			if (pstart[i + 1] - pstart[i] < degree_threshold)
				printf("!!! WA degree in check_core_pruning\n");
			for (ept j = pstart[i]; j < pstart[i + 1]; j++)
				if (edges[j] >= n)
					printf("!!! WA edge in check_core_pruning\n");
		}
	}

	void check_truss_pruning(ui n, ui m, ui triangle_threshold, ept *pstart, ui *edges, char *exists)
	{
		memset(exists, 0, sizeof(char) * n);
		for (ui i = 0; i < n; i++)
		{
			for (ept j = pstart[i]; j < pstart[i + 1]; j++)
				exists[edges[j]] = 1;
			for (ept j = pstart[i]; j < pstart[i + 1]; j++)
			{
				ui v = edges[j], cn = 0;
				for (ept k = pstart[v]; k < pstart[v + 1]; k++)
					if (exists[edges[k]])
						++cn;
				if (cn < triangle_threshold)
					printf("!!! WA triangle in check_truss_pruning\n");
			}
			for (ept j = pstart[i]; j < pstart[i + 1]; j++)
				exists[edges[j]] = 0;
		}
	}
}

#endif /* CTPRUNE_H_ */
