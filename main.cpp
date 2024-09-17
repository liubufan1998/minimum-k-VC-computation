#include "Graph.h"
#include "Utility.h"
#include "Timer.h"
#include "popl.hpp"

using namespace std;
using namespace popl;

void print_usage()
{
	printf("Example usage: ./VCtoB -g path_to_graph -k 10\n");
}

int main(int argc, char *argv[])
{
#ifndef NDEBUG
	printf("**** VCtoB (Debug) build at %s %s ***\n", __TIME__, __DATE__);
	printf("!!! You may want to define NDEBUG in Utility.h to get better performance!\n");
#else
	printf("**** VCtoB (Release) build at %s %s ***\n", __TIME__, __DATE__);
#endif

	OptionParser op("Allowed options");
	auto help_option = op.add<Switch>("h", "help", "\'produce help message\'");
	auto graph_option = op.add<Value<string>>("g", "graph", "\'path to input graph file\'");
	auto k_option = op.add<Value<int>>("k", "k", "\'the value of k for k-VC\'");

	op.parse(argc, argv);

	if (help_option->is_set() || argc <= 1)
	{
		cout << op << endl;
		if (argc <= 1)
		{
			print_usage();
			return 0;
		}
	}

	if (!graph_option->is_set())
	{
		printf("!!! The argument -g is required! Exit !!!\n");
		print_usage();
		return 0;
	}

	if (!k_option->is_set())
	{
		printf("!!! The argument -k is required! Exit !!!\n");
		print_usage();
		return 0;
	}

	Graph *graph = new Graph(graph_option->value().c_str(), k_option->value());
	graph->readBinFile();

	int result = graph->get_min_kVC(); // the variable result records the size of the minimum k-VC. if result == -1 means no k-VC

	delete graph;
	return 0;
}
