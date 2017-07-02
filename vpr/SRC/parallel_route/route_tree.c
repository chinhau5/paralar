#include "route_tree.h"
#include "log.h"
#include "utility.h"

using namespace std;

void route_tree_init(route_tree_t &rt, const GlobalRRGraph *rrg)
{
	rt.rrg = rrg;
	rt.root_rt_node_id = RouteTree::null_vertex();
	rt.num_nodes = 0;
}

void route_tree_clear(route_tree_t &rt)
{
	rt.root_rt_nodes.clear();
	rt.root_rt_node_id = RouteTree::null_vertex();
	rt.num_nodes = 0;
	clear_vertices(rt.graph);
	assert(num_vertices(rt.graph) == 0 && num_edges(rt.graph) == 0);
	rt.rr_node_to_rt_node.clear();
}

bool route_tree_empty(const route_tree_t &rt)
{
	/*int num_rt_nodes = 0;*/
	/*for (const auto &n : route_tree_get_nodes(rt)) {*/
		/*++num_rt_nodes;*/
	/*}*/

	int num_rt_edges = num_edges(rt.graph);
	/*assert((num_rt_nodes > 0 && num_rt_edges > 0)*/
			/*|| (num_rt_nodes == 0 && num_rt_edges == 0));*/
	/*if (rt.root_rt_node_id == -1) {*/
		/*assert(num_rt_nodes == 0);*/
	/*} else {*/
		/*assert(num_rt_nodes > 0);*/
	/*}*/
	/*return rt.root_rt_node_id == -1;*/

	/*assert(num_rt_nodes == rt.num_nodes);*/

	return rt.num_nodes == 0 && num_rt_edges == 0;
}

int route_tree_num_nodes(const route_tree_t &rt)
{
	return rt.num_nodes;
}

boost::iterator_range<boost::filter_iterator<valid_rt_node, route_tree_t::vertex_iterator>>
route_tree_get_nodes(const route_tree_t &rt)
{
	valid_rt_node validator(rt.graph);
	const auto &nodes = get_vertices(rt.graph);

	return boost::make_iterator_range(
			boost::make_filter_iterator<valid_rt_node>(validator, begin(nodes), end(nodes)),
			boost::make_filter_iterator<valid_rt_node>(validator, end(nodes), end(nodes))
			);
	/*return get_vertices(rt.graph) | boost::adaptors::filtered([&rt] (unsigned long rt_node) -> bool { return get_vertex_props(rt.graph, rt_node).valid; }); */
	/*return route_tree_get_nodes_impl<const route_tree_t, route_tree_t::vertex_iterator>(rt);*/
}


boost::iterator_range<route_tree_t::branch_iterator>
route_tree_get_branches(const route_tree_t &rt, RouteTreeNode rt_node)
{
	return get_out_edges(rt.graph, rt_node);
}

RouteTreeNode route_tree_add_rr_node(route_tree_t &rt, GlobalRRNode rr_node)
{
	const auto &iter = rt.rr_node_to_rt_node.find(rr_node);

	RouteTreeNode rt_node;
	rt_node_property_t *rt_node_p;

	if (iter == rt.rr_node_to_rt_node.end()) {
		add_vertex(rt.graph);
		rt_node = num_vertices(rt.graph)-1;
		rt.rr_node_to_rt_node[rr_node] = rt_node;
		rt_node_p = &get_vertex_props(rt.graph, rt_node);
		rt_node_p->rr_node = rr_node;
		/* a sneaky bug when the route tree is cleared and new nodes reuse the old memory */
		/* causing the check assert(!rt_node_p->valid) to fail */
		rt_node_p->valid = false;
	} else {
		rt_node_p = &get_vertex_props(rt.graph, iter->second);
		assert(rt_node_p->rr_node == rr_node);
		if (rt_node_p->valid) {
			rt_node = RouteTree::null_vertex();
		} else {
			rt_node = iter->second;
		}
	}

	if (rt_node != RouteTree::null_vertex()) {
		assert(!rt_node_p->valid);
		rt_node_p->valid = true;
		rt_node_p->pending_rip_up = false;
		rt_node_p->ripped_up = false;
		rt_node_p->rt_edge_to_parent = RouteTree::null_edge();
		/*rt_node_p->branch_point = false;*/
		rt_node_p->num_iterations_fixed = 0;

		rt_node_p->reexpand = false;
		rt_node_p->upstream_R = std::numeric_limits<float>::max();
		rt_node_p->delay = std::numeric_limits<float>::max();

		/*const auto &rr_node_p = get_vertex_props(g, rr_node);*/

		/*segment seg(point(rr_node_p.xlow, rr_node_p.ylow), point(rr_node_p.xhigh, rr_node_p.yhigh));*/

		/*if (rr_node_p.type != IPIN && rr_node_p.type != SINK) {*/
			/*assert(rt.point_tree.count(make_pair(seg, rt_node)) == 0);*/

			/*rt.point_tree.insert(make_pair(seg, rt_node));*/
		/*}*/

		++rt.num_nodes;

		/*assert(rt.num_nodes == rt.point_tree.size());*/
	}

	return rt_node;
}

const RouteTreeEdge &route_tree_add_edge_between_rr_node(route_tree_t &rt, GlobalRRNode rr_node_a, GlobalRRNode rr_node_b)
{
	/*int rt_node_a = route_tree_add_rr_node(rt, rr_node_a);*/
	/*int rt_node_b = route_tree_add_rr_node(rt, rr_node_b);*/
	/*assert(rt_node_a < num_vertices(rt.graph));*/
	/*assert(rt_node_b < num_vertices(rt.graph));*/
	RouteTreeNode rt_node_a = route_tree_get_rt_node(rt, rr_node_a);
	RouteTreeNode rt_node_b = route_tree_get_rt_node(rt, rr_node_b);

	assert(rt_node_a != RouteTree::null_vertex());
	assert(rt_node_b != RouteTree::null_vertex());

	if (has_edge(rt.graph, rt_node_a, rt_node_b)) {
		char buffer_a[256];
		char buffer_b[256];
		sprintf_rr_node(get_vertex_props(rt.graph, rt_node_a).rr_node, buffer_a);
		sprintf_rr_node(get_vertex_props(rt.graph, rt_node_b).rr_node, buffer_b);
		zlog_error(delta_log, "Existing edge between %s and %s\n", buffer_a, buffer_b);
		assert(false);
	}

	const auto &edge = add_edge(rt.graph, rt_node_a, rt_node_b);
	/* a route tree node can only have one driver */
	auto &rt_node_b_p = get_vertex_props(rt.graph, rt_node_b);
	assert(!valid(rt_node_b_p.rt_edge_to_parent));
	rt_node_b_p.rt_edge_to_parent = edge;

	return edge;
}

bool route_tree_has_edge(const route_tree_t &rt, GlobalRRNode a, GlobalRRNode b)
{
	RouteTreeNode rt_node_a = route_tree_get_rt_node(rt, a);
	RouteTreeNode rt_node_b = route_tree_get_rt_node(rt, b);

	assert(rt_node_a != RouteTree::null_vertex());
	assert(rt_node_b != RouteTree::null_vertex());

	return has_edge(rt.graph, rt_node_a, rt_node_b);
}

RouteTreeNode route_tree_get_rt_node(const route_tree_t &rt, GlobalRRNode rr_node)
{
	RouteTreeNode res;

	const auto &iter = rt.rr_node_to_rt_node.find(rr_node);

	if (iter == rt.rr_node_to_rt_node.end()) {
		res = RouteTree::null_vertex();
	} else {
		const auto &v = get_vertex_props(rt.graph, iter->second);
		if (v.valid) {
			res = iter->second;
		} else {
			res = RouteTree::null_vertex();
		}
	}

	return res;
}

void route_tree_set_node_properties(route_tree_t &rt, const RouteTreeNode &rt_node, bool reexpand, float upstream_R, float delay)
{
	auto &rt_node_p = get_vertex_props(rt.graph, rt_node);

	rt_node_p.reexpand = reexpand;
	/*rt_node_p.rr_edge_to_parent = prev_edge;*/
	rt_node_p.upstream_R = upstream_R;
	rt_node_p.delay = delay;
}

void route_tree_add_root(route_tree_t &rt, GlobalRRNode rr_node)
{
	const auto &rt_node = route_tree_get_rt_node(rt, rr_node);
	assert(rt_node != RouteTree::null_vertex());
	assert(find(begin(rt.root_rt_nodes), end(rt.root_rt_nodes), rt_node) == end(rt.root_rt_nodes));

	rt.root_rt_nodes.push_back(rt_node);
}

void route_tree_set_root(route_tree_t &rt, GlobalRRNode rr_node)
{
	const auto &rt_node = route_tree_get_rt_node(rt, rr_node);
	assert(rt_node != RouteTree::null_vertex());
	assert(rt.root_rt_node_id == RouteTree::null_vertex());
	rt.root_rt_node_id = rt_node;
}

void route_tree_remove_edge(route_tree_t &rt, const RouteTreeEdge &rt_edge)
{
	RouteTreeNode from = get_source(rt.graph, rt_edge);
	RouteTreeNode to = get_target(rt.graph, rt_edge);

	const auto &from_p = get_vertex_props(rt.graph, from);
	auto &to_p = get_vertex_props(rt.graph, to);

	char s_from[256];
	char s_to[256];
	sprintf_rr_node(from_p.rr_node, s_from);
	sprintf_rr_node(to_p.rr_node, s_to);
	zlog_level(delta_log, ROUTER_V2, "Removing edge %s -> %s\n", s_from, s_to);

	remove_edge(rt.graph, rt_edge);
	assert(valid(to_p.rt_edge_to_parent));
	to_p.rt_edge_to_parent = RouteTree::null_edge();
}

void route_tree_remove_node(route_tree_t &rt, GlobalRRNode rr_node)
{
	RouteTreeNode rt_node = route_tree_get_rt_node(rt, rr_node);
	auto &rt_node_p = get_vertex_props(rt.graph, rt_node);

	assert(rt_node_p.valid);
	rt_node_p.valid = false;

	/*if (rt_node == rt.root_rt_node_id) {*/
		/*rt.root_rt_node_id = RouteTree::null_vertex();*/
	/*}*/
	const auto &root = find(begin(rt.root_rt_nodes), end(rt.root_rt_nodes), rt_node);
	if (root != end(rt.root_rt_nodes)) {
		rt.root_rt_nodes.erase(root);
	}

	/*const auto &rr_node_p = get_vertex_props(g, rr_node);*/

	/*if (rr_node_p.type != IPIN && rr_node_p.type != SINK) {*/
		/*bool success = rt.point_tree.remove(make_pair(*/
						/*segment(*/
							/*point(rr_node_p.xlow, rr_node_p.ylow),*/
							/*point(rr_node_p.xhigh, rr_node_p.yhigh)*/
							/*),*/
						/*rt_node*/
						/*)*/
					/*);*/

		/*assert(success);*/
	/*}*/

	--rt.num_nodes;
	
	/*assert(rt.num_nodes == rt.point_tree.size());*/
}
