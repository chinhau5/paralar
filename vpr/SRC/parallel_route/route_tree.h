#ifndef ROUTE_TREE_H
#define ROUTE_TREE_H

#include "log.h"
#include "utility.h"
#include "geometry.h"
#include "new_rr_graph.h"
#include <boost/range/iterator_range.hpp>
#include <boost/range/adaptor/filtered.hpp>

typedef struct rt_edge_property_t {
	GlobalRREdge rr_edge;
} rt_edge_property_t;

typedef cache_edge_t<rt_edge_property_t> RouteTreeEdge;

typedef struct rt_node_property_t {
	GlobalRRNode rr_node;
	bool valid;
	bool pending_rip_up;
	bool ripped_up;
	//bool branch_point;
	RouteTreeEdge rt_edge_to_parent;
	int num_iterations_fixed;
	//int saved_num_out_edges; [> only valid for SOURCE <]

	//int owner;

	/* public properties */
	bool reexpand;
	//RREdge rr_edge_to_parent;
	float upstream_R;	
	float delay;
	float downstream_C;
	/*float upstream_R_from_route_state;*/
} rt_node_property_t;

typedef cache_graph_t<rt_node_property_t, rt_edge_property_t> RouteTree;
typedef int RouteTreeNode;

typedef struct route_tree_t {
	typedef std::pair<segment, int> rtree_value;

	const GlobalRRGraph *rrg;
	RouteTree graph;
	std::vector<int> root_rt_nodes;
	int root_rt_node_id;
	int num_nodes;
	std::map<int, int> rr_node_to_rt_node;
	//map<int, vector<int>> sink_rr_node_to_path;
	//map<int, RouteTreeNode> path_branch_point;
	//bgi::rtree<rtree_value, bgi::rstar<64>> point_tree;
	//box scheduler_bounding_box;
	//std::map<int, std::vector<int>> sink_edges;

	//template<typename Graph, typename Value, typename Base>
	//struct iterator : public std::iterator<
					  //typename std::forward_iterator_tag,
					  //Value,
					  //ptrdiff_t,
					  //Value,
					  //Value 
					  //> { 

		//const Graph &g;
		//Base c;
		//Base e;

		//iterator(const Graph &g, const Base &c, const Base &e) 
			//: g(g), c(c), e(e)
		   //{
		//}

		//iterator &operator++()
		//{
			//do {
				//if (c != e) {
					//++c;
				//}
			//} while (c != e && get_vertex(g, c).properties.valid);
		
			//return *this;
		//}

		//typename iterator::reference operator*() const
		//{
			//return c;
		//}

		//bool operator==(const iterator &other) const
		//{
			//return other.c == c;
		//}

		//bool operator!=(const iterator &other) const
		//{
			//return other.c != c;
		//}
	//};

	//typedef iterator<RouteTree, int, int> vertex_iterator;
	typedef typename RouteTree::vertex_iterator vertex_iterator;
	//typedef iterator<const RouteTreeNode, typename RouteTree::vertex_iterator> vertex_const_iterator;
	typedef typename RouteTree::out_edge_iterator branch_iterator;
	//typedef typename RouteTree::out_edges_const_iterator branch_const_iterator;
} route_tree_t;

void route_tree_init(route_tree_t &rt, const GlobalRRGraph *rrg);

void route_tree_clear(route_tree_t &rt);

bool route_tree_empty(const route_tree_t &rt);

int route_tree_num_nodes(const route_tree_t &rt);

RouteTreeNode route_tree_add_rr_node(route_tree_t &rt, GlobalRRNode rr_node);

const RouteTreeEdge &route_tree_add_edge_between_rr_node(route_tree_t &rt, GlobalRRNode rr_node_a, GlobalRRNode rr_node_b);

void route_tree_remove_edge(route_tree_t &rt, const RouteTreeEdge &rt_edge);

void route_tree_remove_node(route_tree_t &rt, GlobalRRNode rr_node);

bool route_tree_has_edge(const route_tree_t &rt, GlobalRRNode a, GlobalRRNode b);

RouteTreeNode route_tree_get_rt_node(const route_tree_t &rt, GlobalRRNode rr_node);

struct valid_rt_node {
	const RouteTree &g;
	valid_rt_node(const RouteTree &g) :
		g(g) {}

	bool operator()(unsigned long rt_node) const
	{
		return get_vertex_props(g, rt_node).valid; 
	}
};

boost::iterator_range<boost::filter_iterator<valid_rt_node, route_tree_t::vertex_iterator>>
route_tree_get_nodes(const route_tree_t &rt);

boost::iterator_range<route_tree_t::branch_iterator>
route_tree_get_branches(const route_tree_t &rt, RouteTreeNode rt_node);

void route_tree_set_root(route_tree_t &rt, GlobalRRNode rr_node);

void route_tree_add_root(route_tree_t &rt, GlobalRRNode rr_node);

void route_tree_set_node_properties(route_tree_t &rt, const RouteTreeNode &rt_node, bool reexpand, float upstream_R, float delay);

//template<typename Congestion>
//bool route_tree_is_congested_internal(route_tree_t &rt, RouteTreeNode rt_node, const RRGraph &g, const Congestion *congestion)
//{
	//bool is_congested = route_tree_node_is_congested(rt, rt_node, g, congestion);

	//if (!is_congested) {
		//auto bis = route_tree_get_branches(rt, rt_node);

		//for (auto bi = std::begin(bis); !is_congested && bi != std::end(bis); ++bi) {
			//const auto &child = get_target(rt.graph, *bi);

			//is_congested = route_tree_is_congested_internal(rt, child, g, congestion);
		//}
	//}

	//return is_congested;
//}

//template<typename Congestion>
//bool route_tree_is_congested(route_tree_t &rt, const RRGraph &g, const Congestion *congestion)
//{
	//int num_marked = 0;

	//assert(rt.root_rt_nodes.size() == 1 || rt.root_rt_nodes.size() == 0);

	//bool is_congested = false;
	//for (int i = 0; i < rt.root_rt_nodes.size() && !is_congested; ++i) {
		//is_congested = route_tree_is_congested_internal(rt, rt.root_rt_nodes[i], g, congestion);
	//}

	//return is_congested;
//}

#endif
