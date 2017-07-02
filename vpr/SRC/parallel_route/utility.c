#include <assert.h>
#include "vpr_types.h"
#include "utility.h"

static const char *rr_types[] =  {
	"SOURCE", "SINK", "IPIN", "OPIN", "CHANX", "CHANY", "INTRA_CLUSTER_EDGE"
};

/*#define PRINT_RR_NODE*/

static const GlobalRRGraph *rr_graph;

void init_sprintf_rr_node(const GlobalRRGraph *_rr_graph)
{
	rr_graph = _rr_graph;
}

void sprintf_rr_node_impl(int inode, char *buffer)
{
	const auto &props = get_vertex_props(*rr_graph, inode);

	sprintf(buffer, "%d %s (%d,%d)(%d,%d) ", inode, rr_types[props.type], props.xlow, props.ylow, props.xhigh, props.yhigh);
}
