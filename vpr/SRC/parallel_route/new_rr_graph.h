#ifndef NEW_RR_GRAPH_H
#define NEW_RR_GRAPH_H

#include "vpr_types.h"
#include "cache_graph.h"

typedef struct rr_node_property_t {
	t_rr_type type;
	bool inc_direction;
	int xlow;
	int ylow;
	int xhigh;
	int yhigh;
	//int pin_height;
	int real_xlow;
	int real_ylow;
	int real_xhigh;
	int real_yhigh;
	float R;
	float C;
	int cost_index;
	int capacity;
	//int occ;
	//int recalc_occ;
	//float pres_cost;
	//float acc_cost;
	//tbb::spin_mutex *lock;
} rr_node_property_t;

//typedef struct rr_node_property_t {
	//t_rr_type type;
	//bool inc_direction;
	//int xlow;
	//int ylow;
	//int xhigh;
	//int yhigh;
	//int real_xlow;
	//int real_ylow;
	//int real_xhigh;
	//int real_yhigh;
	//float R;
	//float C;
	//int cost_index;
	//tbb::spin_mutex *lock;
	//int capacity;
	//int occ;
	//std::vector<int> users;
	//int recalc_occ;
	//float pres_cost;
	//float acc_cost;
//} rr_node_property_t;

typedef struct rr_edge_property_t {
	int id;
	int switch_index;
	//bool buffered;
	//float switch_delay;
	//float R;
} rr_edge_property_t;

typedef cache_graph_t<rr_node_property_t, rr_edge_property_t> RRGraph;
typedef int RRNode;
typedef cache_edge_t<rr_edge_property_t> RREdge;

struct global_rr_node_property_t {
	t_rr_type type;
	int xlow;
	int ylow;
	int xhigh;
	int yhigh;
};

struct global_rr_edge_property_t {
	bool ignore;
	int util;
	int capacity;
	float mult;
};

typedef cache_graph_t<global_rr_node_property_t, global_rr_edge_property_t> GlobalRRGraph;
typedef int GlobalRRNode;
typedef cache_edge_t<global_rr_edge_property_t> GlobalRREdge;

#endif
