#include <utility>
#include <queue>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include "util.h"
#include "vpr_types.h"
#include "globals.h"
#include "route_export.h"
#include "route_common.h"
#include "route_tree_timing.h"
#include "route_timing.h"
#include "heapsort.h"
#include "path_delay.h"
#include "net_delay.h"
#include "stats.h"
#include "ReadOptions.h"

#include "new_rr_graph.h"
#include "geometry.h"
#include "log.h"
#include "dijkstra.h"
#include "route_tree.h"

using namespace std;

struct source_t {
	struct net_t *net;
	int rr_node;
	int x;
	int y;
};

struct sink_t {
	struct net_t *net;
	int id;
	float criticality_fac;
	int rr_node;
	int x;
	int y;
	box bounding_box;
};

struct net_t {
	int vpr_id;
	int local_id;
	source_t source;
	std::vector<sink_t> sinks;
	box bounding_box;
};

static bool inside_bb(const global_rr_node_property_t &node, const box &bb)
{
	//int xlow, xhigh, ylow, yhigh;

	//switch (node.type) {
		//case CHANX:
		//case CHANY:
			//if (node.inc_direction) {
				//xlow = node.xlow;
				//ylow = node.ylow;
			//} else {
				//xlow = node.xhigh;
				//ylow = node.yhigh;
			//}
			//xhigh = xlow;
			//yhigh = ylow;
			//break;

		//case IPIN:
		//case OPIN:
			//xlow = node.xlow;
			//ylow = node.ylow;
			//xhigh = node.xhigh;
			//yhigh = node.ylow + node.pin_height;
			//assert(xlow == xhigh);
			//assert(ylow <= yhigh);
			//break;

		//case SOURCE:
		//case SINK:
			//xlow = node.xlow;
			//ylow = node.ylow;
			//xhigh = node.xhigh;
			//yhigh = node.yhigh;
			//break;

		//default:
			//assert(false);
			//break;
	//}

	bool inside;

	if (node.xhigh < bg::get<bg::min_corner, 0>(bb)
			|| node.xlow > bg::get<bg::max_corner, 0>(bb)
			|| node.yhigh < bg::get<bg::min_corner, 1>(bb)
			|| node.ylow > bg::get<bg::max_corner, 1>(bb)) {
		inside = false;
	} else {
		inside = true;
	}

	return inside;
}

static bool inside_bbs(const global_rr_node_property_t &node, const vector<box> &boxes)
{
	bool inside = false;
	for (int i = 0; i < boxes.size() && !inside; ++i) {
		inside = inside_bb(node, boxes[i]);
	}
	return inside;
}


tuple<int, int, int, int, int> get_tuple(int inode)
{
	extern s_rr_node *rr_node;

	int type, real_xlow, real_ylow, real_xhigh, real_yhigh;

	//extern struct s_grid_tile **grid;
	//auto tile = &grid[rr_node[inode].xlow][rr_node[inode].ylow];
	//int pin_height;

	switch (rr_node[inode].type) {
		case CHANX:
		case CHANY:
			type = rr_node[inode].type;
			if (rr_node[inode].direction == INC_DIRECTION) {
				real_xlow = rr_node[inode].xlow;
				real_ylow = rr_node[inode].ylow;
				real_xhigh = rr_node[inode].xhigh;
				real_yhigh = rr_node[inode].yhigh;
			} else {
				real_xlow = rr_node[inode].xhigh;
				real_ylow = rr_node[inode].yhigh;
				real_xhigh = rr_node[inode].xlow;
				real_yhigh = rr_node[inode].ylow;
			}
			break;

		case IPIN:
		case OPIN:
			type = -1;
			real_xlow = -1;
			real_ylow = -1;
			real_xhigh = -1;
			real_yhigh = -1;
			//type = rr_node[inode].type;
			//num_nodes = 1;
			//real_xlow = rr_node[inode].xlow;
			//pin_height = tile->type->pin_height[rr_node[inode].ptc_num];
			//real_ylow = rr_node[inode].ylow + pin_height;
			//assert(rr_node[inode].xlow == rr_node[inode].xhigh);

			//real_xhigh = real_xlow;
			//real_yhigh = real_ylow;

			break;

		case SOURCE:
		case SINK:
			type = rr_node[inode].type;
			real_xlow = rr_node[inode].xlow;
			real_ylow = rr_node[inode].ylow;
			real_xhigh = rr_node[inode].xhigh;
			real_yhigh = rr_node[inode].yhigh;

			break;

		default:
			assert(false);
			break;
	}

	return make_tuple(type, real_xlow, real_ylow, real_xhigh, real_yhigh);
}

template<typename Callback>
void bfs(int root_inode, const Callback &callback)
{
	extern s_rr_node *rr_node;

	struct qitem {
		int inode;
		int level;
	};

	queue<qitem> q;

	//int level;
	//switch () {
		//case CHANX:
		//case CHANY:
			//level = 
			//break;

		//case SOURCE:
		//case SINK:
			//break;
	//}
	q.push({ root_inode, 2 });

	while (!q.empty()) {
		qitem item = q.front();
		q.pop();

		int type = rr_node[item.inode].type;

		if (item.level == 0 || (item.level == 1 && (type == CHANX || type == CHANY))) {
			assert(type == CHANX || type == CHANY || type == SINK);
			callback(item.inode);
		} else {
			for (int i = 0; i < rr_node[item.inode].num_edges; ++i) {
				int v = rr_node[item.inode].edges[i];
				q.push({ v, item.level-1 });
			}
		}
	}
}

struct extra_route_state_t {
	int same;
	int ortho;
	bool existing;
	float upstream_R;
	float delay;
};

struct dijkstra_stats_t {
	unsigned long num_heap_pops;
	unsigned long num_heap_pushes;
	unsigned long num_neighbor_visits;
};

struct congestion_t {
	//tbb::spin_mutex lock;
	int occ;
	float pres_cost;
	float acc_cost;
	int recalc_occ;
};

class SinkRouter {
	private:
		const GlobalRRGraph &_g;

		float _astar_fac;

		float *_known_distance;
		float *_distance;
		GlobalRREdge *_prev_edge;
		extra_route_state_t *_state;

		vector<int> _modified_nodes;
		vector<bool> _modified_node_added;

		GlobalRRNode _existing_opin;

		const sink_t *_current_sink;
		const vector<box> *m_bounding_boxes;

		dijkstra_stats_t _stats;

	private:
		void popped_node(const heap_node_t<GlobalRREdge, extra_route_state_t> &node)
		{
			char buffer[256];
			sprintf_rr_node(node.node, buffer);
			zlog_level(delta_log, ROUTER_V3, "Current: %s Existing: %d [kd: %g %a okd: %g %a] [d: %g %a od: %g %a] prev: %d\n",
					buffer, node.extra.existing,
					node.known_distance, node.known_distance,
					_known_distance[node.node], _known_distance[node.node],
					node.distance, node.distance,
					_distance[node.node], _distance[node.node],
					get_source(_g, node.prev_edge));

			++_stats.num_heap_pops;
		}

		void relax_node(const heap_node_t<GlobalRREdge, extra_route_state_t> &node)
		{
			//RouteTreeNode rt_node = route_tree_get_rt_node(*_current_rt, v);

			//if (rt_node != RouteTree::null_vertex()) {
				//const auto &rt_node_p = get_vertex_props(_current_rt->graph, rt_node);
				//if (!valid(_prev_edge[v])) {
					//assert(!valid(rt_node_p.rt_edge_to_parent));
				//} else {
					//int rr = get_vertex_props(_current_rt->graph, get_source(_current_rt->graph, rt_node_p.rt_edge_to_parent)).rr_node;
					//assert(rr == get_source(_g, _prev_edge[v]));
				//}
			//}

			//int v = node.node;
			//const auto &v_p = get_vertex_props(_g, v);

			//if (valid(node.prev_edge)) {
				//assert(v == get_target(_g, node.prev_edge));

				//int u = get_source(_g, node.prev_edge);

				//float u_delay;
				//float u_upstream_R;
				//if (_state[u].upstream_R != std::numeric_limits<float>::max()) {
					//assert(_state[u].delay != std::numeric_limits<float>::max());
					//u_upstream_R = _state[u].upstream_R;
					//u_delay = _state[u].delay;
				//} else {
					//auto rt_node = route_tree_get_rt_node(*_current_rt, u);
					//assert(rt_node != RouteTree::null_vertex());
					//const auto &u_rt_node_p = get_vertex_props(_current_rt->graph, rt_node);
					//u_upstream_R = u_rt_node_p.upstream_R;
					//u_delay = u_rt_node_p.delay;
				//}

				//const auto &e_p = get_edge_props(_g, node.prev_edge);

				//extern struct s_switch_inf *switch_inf;
				//const struct s_switch_inf *sw = &switch_inf[e_p.switch_index];

				//_state[v].upstream_R = sw->R + v_p.R;
				//if (!sw->buffered)  {
					//_state[v].upstream_R += u_upstream_R;
				//} 

				//float delay;
				//if (sw->buffered) {
					//delay = sw->Tdel + v_p.C * (sw->R + 0.5 * v_p.R);
				//} else {
					//delay = sw->Tdel + v_p.C * (u_upstream_R + sw->R + 0.5 * v_p.R);
				//}
				//_state[v].delay = u_delay + delay;
			//} else {
				//_state[v].upstream_R = v_p.R;
				//_state[v].delay = v_p.C * 0.5 * v_p.R;
			//}
			_state[node.node] = node.extra;

			if (!_modified_node_added[node.node]) {
				_modified_nodes.push_back(node.node);
				_modified_node_added[node.node] = true;
			}
		}

		bool expand_node(int node)
		{
			++_stats.num_neighbor_visits;

			return true;
			
			//const auto &prop = get_vertex_props(_g, node);

			//char buffer[256];
			//sprintf_rr_node(node, buffer);

			//zlog_level(delta_log, ROUTER_V3, "\tChecking whether to expand %s ", buffer);

			//if (!inside_bbs(prop, *m_bounding_boxes)) {
				//zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
				//return false;
			//}

			//zlog_level(delta_log, ROUTER_V3, " OK\n");

			//return true;
		}

		void get_edge_weight(const GlobalRREdge &e, float &known_distance, float &distance, extra_route_state_t &extra)
		{
			int u = get_source(_g, e);
			int v = get_target(_g, e);

			const auto &v_p = get_vertex_props(_g, v);
			const auto &e_p = get_edge_props(_g, e);

			assert(_state[u].upstream_R != std::numeric_limits<float>::max());

			float delay = 1;
			assert(_state[u].delay != std::numeric_limits<float>::max());
			extra.delay = _state[u].delay + delay; 

			extra.existing = false;

			known_distance = delay + get_edge_props(_g, e).mult;

			distance = known_distance;

			unsigned int xe, xk;
			//memcpy(&xe, &expected_cost, sizeof (xe));
			//memcpy(&xk, &known_distance, sizeof (xk));

			zlog_level(delta_log, ROUTER_V3, "\t%d -> %d delay %g known %g %a %X actual %g %a %X\n", 
					u, v, delay, known_distance, known_distance, xk, distance , distance, xe);
			//zlog_level(delta_log, ROUTER_V3, "\t[u: upstream %g] [edge: d %g R %g] [v: R %g C %g]\n",
					//_state[u].upstream_R, sw->Tdel, sw->R, v_p.R, v_p.C);
		}

		void push_node(const heap_node_t<GlobalRREdge, extra_route_state_t> &node)
		{
			zlog_level(delta_log, ROUTER_V3, "\tPushing neighbor %d [kd %g %a d %g %a] to heap\n",
					node.node,
					node.known_distance, node.known_distance,
					node.distance, node.distance);

			++_stats.num_heap_pushes;
		}

		void backtrack(int sink_node, route_tree_t &rt, vector<GlobalRRNode> &added_rr_nodes)
		{
			GlobalRRNode rr_node = sink_node;
			GlobalRRNode child_rr_node = GlobalRRGraph::null_vertex();
			GlobalRREdge edge = GlobalRRGraph::null_edge();

			bool stop = false;
			while (!stop) {
				const auto &rr_node_p = get_vertex_props(_g, rr_node);

				char buffer[256];
				sprintf_rr_node(rr_node, buffer);
				zlog_level(delta_log, ROUTER_V3, "Adding %s to route tree\n", buffer);

				RouteTreeNode rt_node = route_tree_add_rr_node(rt, rr_node);
				if (rt_node != RouteTree::null_vertex()) {
					assert(inside_bbs(rr_node_p, *m_bounding_boxes));

					assert(_state[rr_node].upstream_R != std::numeric_limits<float>::max() && _state[rr_node].delay != std::numeric_limits<float>::max());
					route_tree_set_node_properties(rt, rt_node, rr_node_p.type != IPIN && rr_node_p.type != SINK, _state[rr_node].upstream_R, _state[rr_node].delay);

					if (rr_node_p.type == SOURCE) {
						route_tree_add_root(rt, rr_node);
					}

					added_rr_nodes.push_back(rr_node);
				} else {
					stop = true;
				}

				if (rr_node_p.type == OPIN && _existing_opin == GlobalRRGraph::null_vertex()) {
					_existing_opin = rr_node;
				}

				if (child_rr_node != GlobalRRGraph::null_vertex()) {
					char buffer2[256];
					sprintf_rr_node(child_rr_node, buffer2);

					zlog_level(delta_log, ROUTER_V3, "Adding edge %s -> %s\n", buffer, buffer2);

					const RouteTreeEdge &rt_edge = route_tree_add_edge_between_rr_node(rt, rr_node, child_rr_node); 
					auto &rt_edge_props = get_edge_props(rt.graph, rt_edge);
					assert(valid(edge));
					assert(get_source(_g, edge) == rr_node && get_target(_g, edge) == child_rr_node);
					rt_edge_props.rr_edge = edge;
				}

				zlog_level(delta_log, ROUTER_V3, "\n");

				child_rr_node = rr_node;
				edge = _prev_edge[rr_node];
				if (valid(edge)) {
					rr_node = get_source(_g, edge);
				} else {
					stop = true;
				}
			}
		}

		void reset_stats()
		{
			_stats.num_heap_pops = 0;
			_stats.num_heap_pushes = 0;
			_stats.num_neighbor_visits = 0;
		}

		//template<typename Graph, typename Edge, typename EdgeWeightFunc, typename Callbacks>
		//friend void delta_stepping(const Graph &g, const vector<heap_node_t<Edge>> &sources, int sink, float delta, float *known_distance, float *distance, Edge *prev_edge, const EdgeWeightFunc &edge_weight, Callbacks &callbacks);

		template<typename Graph, typename Edge, typename EdgeWeightFunc, typename ExpandCheckFunc, typename Callbacks, typename Extra>
		friend void dijkstra(const Graph &g, priority_queue<heap_node_t<Edge, Extra>> &sources, int sink, float *known_distance, float *distance, Edge *prev_edge, const ExpandCheckFunc &expand_node, const EdgeWeightFunc &edge_weight, Callbacks &callbacks);

		//template<typename Edge, typename Callbacks>
		//friend void relax(Buckets &buckets, float delta, vector<bool> &in_bucket, const vector<bool> &vertex_deleted,
					//float *known_distance, float *distance, Edge *predecessor,
					//int v, float new_known_distance, float new_distance, const Edge &edge,
					//Callbacks &callbacks);

	public:
		SinkRouter(const GlobalRRGraph &g)
			: _g(g), _modified_node_added(num_vertices(g), false) 
		{
			_known_distance = new float[num_vertices(g)];
			_distance = new float[num_vertices(g)];
			_prev_edge = new GlobalRREdge[num_vertices(g)];
			_state = new extra_route_state_t[num_vertices(g)];

			for (int i = 0; i < num_vertices(g); ++i) {
				_known_distance[i] = std::numeric_limits<float>::max();
				_distance[i] = std::numeric_limits<float>::max();
				_prev_edge[i] = GlobalRRGraph::null_edge();
				_state[i].delay = std::numeric_limits<float>::max();
				_state[i].upstream_R = std::numeric_limits<float>::max();
			}

			reset_stats();
		}

		void route(const source_t *source, const sink_t *sink, const vector<box> &bounding_boxes, float astar_fac, route_tree_t &rt, vector<GlobalRRNode> &added_rr_nodes, float &delay, dijkstra_stats_t *stats)
		{
			_current_sink = sink;
			m_bounding_boxes = &bounding_boxes;
			_astar_fac = astar_fac;

			/* TODO: we should not reset the exisiting OPIN here */
			_existing_opin = GlobalRRGraph::null_vertex();

			reset_stats();

			char buffer[256];
			sprintf_rr_node(sink->rr_node, buffer);
			//zlog_level(delta_log, ROUTER_V3, "Current sink: %s BB %d->%d, %d->%d\n", buffer, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax);
			const auto &sink_p = get_vertex_props(_g, sink->rr_node);
			//assert(sink_p.real_xlow == sink_p.real_xhigh && sink_p.real_ylow == sink_p.real_yhigh);
			zlog_level(delta_log, ROUTER_V3, "Current sink: %s Source (%d, %d) BBs ",
					buffer, source ? source->x : -1, source ? source->y : -1);
#ifdef PRINT_BBS 
			for (int i = 0; i < bounding_boxes.size(); ++i) {
				zlog_level(delta_log, ROUTER_V3, "%d %d, %d %d -- ",
						bg::get<bg::min_corner, 0>(bounding_boxes[i]), bg::get<bg::max_corner, 0>(bounding_boxes[i]),
						bg::get<bg::min_corner, 1>(bounding_boxes[i]), bg::get<bg::max_corner, 1>(bounding_boxes[i]));
			}
#endif
			zlog_level(delta_log, ROUTER_V3, "\n");

			assert(source);

			std::priority_queue<heap_node_t<GlobalRREdge, extra_route_state_t>> sources;

			float kd = 0;
			float d = kd; 

			sources.push({ source->rr_node, kd, d, GlobalRRGraph::null_edge(), { 0, 0, false, 0, 0 } });

			sprintf_rr_node(source->rr_node, buffer);
			zlog_level(delta_log, ROUTER_V3, "Using %s as initial heap\n", buffer);

			//delta_stepping(_g, sources, sink->rr_node, delta, _known_distance, _distance, _prev_edge, [this] (const RREdge &e) -> pair<float, float> { return get_edge_weight(e); }, *this);

			dijkstra(_g, sources, sink->rr_node, _known_distance, _distance, _prev_edge, [this] (int rr_node) -> bool { return expand_node(rr_node); }, [this] (const GlobalRREdge &e, float &known_distance, float &distance, extra_route_state_t &extra) -> void { return get_edge_weight(e, known_distance, distance, extra); }, *this);

			backtrack(sink->rr_node, rt, added_rr_nodes);

			assert(get_vertex_props(rt.graph, route_tree_get_rt_node(rt, sink->rr_node)).delay == _state[sink->rr_node].delay);
			//net_timing.delay[sink->id+1] = _state[sink->rr_node].delay;
			delay = _state[sink->rr_node].delay;

			if (stats) {
				*stats = _stats;
			}

			for (const auto &n : _modified_nodes) {
				_known_distance[n] = std::numeric_limits<float>::max();
				_distance[n] = std::numeric_limits<float>::max();
				_prev_edge[n] = GlobalRRGraph::null_edge();
				_state[n].delay = std::numeric_limits<float>::max();
				_state[n].upstream_R = std::numeric_limits<float>::max();

				assert(_modified_node_added[n]);
				_modified_node_added[n] = false;
			}

			_modified_nodes.clear();

			//for (const auto &n : get_vertices(_g)) {
			//assert(_known_distance[n] == std::numeric_limits<float>::max() &&
			//_distance[n] == std::numeric_limits<float>::max() &&
			//_prev_edge[n] == RRGraph::null_edge() &&
			//_state[n].delay == std::numeric_limits<float>::max() &&
			//_state[n].upstream_R == std::numeric_limits<float>::max());
			//}
		}
};

void build_global_graph(GlobalRRGraph &g, map<std::tuple<int, int, int, int, int>, vector<int>> &lookup)
{
	extern s_rr_node *rr_node;
	extern int num_rr_nodes;
	extern t_rr_indexed_data *rr_indexed_data;
	extern struct s_switch_inf *switch_inf;

	for (int i = 0; i < num_rr_nodes; ++i) {
		auto t = get_tuple(i);
		int type = get<0>(t);

		if (type == -1) {
			continue;
		}

		int num_nodes;

		switch (type) {
			case CHANX:
			case CHANY:
				num_nodes = 2;
				break;

			case SOURCE:
			case SINK:
				num_nodes = 1;
				break;

			default:
				assert(false);
				break;
		}

		auto iter = lookup.find(t);
		if (iter == lookup.end()) {
			int base = num_vertices(g);

			add_vertex(g, num_nodes);

			for (int j = 0; j < num_nodes; ++j) {
				auto &v_props = get_vertex_props(g, base+j);

				v_props.type = (t_rr_type)type;
				v_props.xlow = get<1>(t);
				v_props.ylow = get<2>(t);
				v_props.xhigh = get<3>(t);
				v_props.yhigh = get<4>(t);
			}

			vector<int> nodes;
			for (int j = 0; j < num_nodes; ++j) {
				nodes.push_back(base+j);
				if (j < num_nodes-1) {
					const auto &e = add_edge(g, base+j, base+j+1);
					auto &e_props = get_edge_props(g, e);
					e_props.ignore = false;
					e_props.capacity = 0;
				}
			}

			lookup.insert(make_pair(t, nodes));
		} else {
			if (num_nodes == 2) {
				const auto &e = get_edge(g, iter->second[0], iter->second[1]);
				auto &e_props = get_edge_props(g, e);
				++e_props.capacity;
			}

			assert(iter->second.size() == num_nodes);
		}
	}

	for (int i = 0; i < num_rr_nodes; ++i) {
		auto from_t = get_tuple(i);
		auto from_iter = lookup.find(from_t);

		if (from_iter == lookup.end()) {
			assert(rr_node[i].type == IPIN || rr_node[i].type == OPIN);
			continue;
		}

		int u; 
		if (from_iter->second.size() == 2) {
			assert(rr_node[i].type == CHANX || rr_node[i].type == CHANY);
			u = from_iter->second[1];
		} else {
			assert(from_iter->second.size() == 1);
			assert(rr_node[i].type == SOURCE || rr_node[i].type == SINK);
			u = from_iter->second[0];
		}

		bfs(i, [&] (int inode) -> void
				{
					assert(rr_node[inode].type == CHANX || rr_node[inode].type == CHANY ||
							rr_node[inode].type == SOURCE || rr_node[inode].type == SINK);

					auto to_t = get_tuple(inode);
					auto to_iter = lookup.find(to_t);

					assert(to_iter != lookup.end());

					int v = to_iter->second[0];

					if (!has_edge(g, u, v)) {
						const auto &e = add_edge(g, u, v);
						auto &e_props = get_edge_props(g, e);
						e_props.ignore = true;
						e_props.capacity = -1;
					}
				});

		//for (int j = 0; j < rr_node[i].num_edges; ++j) {
			//int to = rr_node[i].edges[j];

			//auto to_t = get_tuple(to);
			//auto to_iter = lookup.find(to_t);

			//if (to_iter == lookup.end()) {
				//assert(rr_node[to].type == IPIN || rr_node[to].type == OPIN);
				//continue;
			//}
			
			//int v;
			//if (to_iter->second.size() == 2) {
				//v = to_iter->second[0];
			//} else {
				//assert(to_iter->second.size() == 1);
				//v = to_iter->second[0];
			//}

			//if (!has_edge(g, u, v)) {
				//add_edge(g, u, v);
			//}
		//}
	}

	using namespace boost::accumulators;
	accumulator_set<int, stats<tag::mean, tag::max, tag::min, tag::variance>> acc;

	for (const auto &e : get_edges(g)) {
		auto &edge_props = get_edge_props(g, e);
		if (!edge_props.ignore) {
			acc(edge_props.capacity);
		}
	}

	printf("global graph has %d vertices %d edges min cap %d max cap %d mean cap %g\n", num_vertices(g), num_edges(g), boost::accumulators::min(acc), boost::accumulators::max(acc), boost::accumulators::mean(acc));
}

void init_nets(vector<net_t> &nets, vector<net_t> &global_nets, const map<std::tuple<int, int, int, int, int>, vector<int>> &lookup, int bb_factor, bool large_bb)
{
	extern s_rr_node *rr_node;
	extern struct s_net *clb_net;
	extern int num_nets;
	extern int **net_rr_terminals; /* [0..num_nets-1][0..num_pins-1] */
	//extern struct s_rr_node *rr_node;
	extern struct s_bb *route_bb;
	extern struct s_block *block;

	int local_id = 0;
	for (int i = 0; i < num_nets; ++i) {
		net_t net;

		net.vpr_id = i;

		int d_rr_node = net_rr_terminals[i][0];

		auto t = get_tuple(d_rr_node);
		auto iter = lookup.find(t);
		assert(iter != lookup.end());
		assert(get<0>(iter->first) == SOURCE);
		assert(iter->second.size() == 1);

		net.source.rr_node = iter->second[0];
		int b = clb_net[i].node_block[0];
		int p = clb_net[i].node_block_pin[0];
		net.source.x = block[b].x;
		net.source.y = block[b].y + block[b].type->pin_height[p];
		//net.source.xhigh = block[b].x;
		//net.source.yhigh = block[b].y + block[b].type->height - 1;
		//assert(net.source.x == rr_node[net.source.rr_node].xlow
				//&& net.source.xhigh == rr_node[net.source.rr_node].xhigh
				//&& net.source.y == rr_node[net.source.rr_node].ylow
				//&& net.source.yhigh == rr_node[net.source.rr_node].yhigh);

		char buffer[256];
		buffer[0] = 0;
		//sprintf_rr_node(net.source.rr_node, buffer);
		//zlog_debug(net_log, "Net %d source %s pin height %d\n", net.vpr_id, buffer, block[b].type->pin_height[p]);

		//net.current_source = net.source;
		/*net.previous_source.rr_node = -1;*/
		/*net.previous_source.x = -1;*/
		/*net.previous_source.y = -1;*/
		/*net.previous_source_valid = false;*/
		 
		for (int j = 1; j <= clb_net[i].num_sinks; j++) {
			sink_t sink;

			d_rr_node = net_rr_terminals[i][j];

			t = get_tuple(d_rr_node);
			iter = lookup.find(t);

			assert(iter != lookup.end());
			assert(get<0>(iter->first) == SINK);
			assert(iter->second.size() == 1);

			sink.rr_node = iter->second[0];
			sink.id = j-1;
			int b = clb_net[i].node_block[j];
			int p = clb_net[i].node_block_pin[j];
			sink.x = block[b].x;
			sink.y = block[b].y + block[b].type->pin_height[p];
			sink.criticality_fac = std::numeric_limits<float>::max();

			//sprintf_rr_node(sink.rr_node, buffer);
			//zlog_debug(net_log, "\tNet %d sink %d (%s) pin height %d, bounding box %d-%d %d-%d\n",
					//net.vpr_id, sink.id, buffer, block[b].type->pin_height[p],
					//sink.current_bounding_box.xmin, sink.current_bounding_box.xmax, sink.current_bounding_box.ymin, sink.current_bounding_box.ymax);
			//sink.congested_iterations = 0;

			net.sinks.push_back(sink);
		}

		assert(net.sinks.size() > 0);

		//if (large_bb) {
			//for (auto &sink : net.sinks) {
				//sink.current_bounding_box = bb;
				//bg::assign_values(sink.bounding_box, bb.xmin, bb.ymin,
						//bb.xmax, bb.ymax);
			//}
		//}

		/*net.box.xmin = route_bb[i].xmin;*/
		/*net.box.ymin = route_bb[i].ymin;*/
		/*net.box.xmax = route_bb[i].xmax;*/
		/*net.box.ymax = route_bb[i].ymax;*/

		if (clb_net[i].is_global) {
			net.local_id = -1;

			global_nets.push_back(net);
		} else {
			net.local_id = local_id;
			++local_id;

			nets.push_back(net);
		}
	}

	/* update pointers */
	for (auto &net : nets) {
		net.source.net = &net;
		for (auto &sink : net.sinks) {
			sink.net = &net;
		}
	}

	printf("Num nets %d num global nets %d\n", nets.size(), global_nets.size());
	/*int num_local_nets = local_id;*/
	/*for (auto &net : nets) {*/
		/*net.num_local_nets = num_local_nets;*/
		/*net.overlapping_nets = new bool[num_local_nets];*/
		/*net.non_overlapping_nets = new bool[num_local_nets];*/
	/*}*/
}

void check_route_tree_internal(const route_tree_t &rt, RouteTreeNode rt_node, GlobalRRGraph &g, vector<int> &visited_sinks, vector<int> &visited_nodes)
{
	const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);

	assert(rt_node_p.valid);

	int rr_node = rt_node_p.rr_node;
	auto &rr_node_p = get_vertex_props(g, rr_node);
	if (rr_node_p.type == SINK) {
		visited_sinks.push_back(rr_node);
		char buffer[256];
		sprintf_rr_node(rr_node, buffer);
		/*zlog_level(delta_log, ROUTER_V2, "route_tree_check: %s\n", buffer);*/
	}

	visited_nodes.push_back(rr_node);

	for (const auto &branch : route_tree_get_branches(rt, rt_node)) {
		/*for_all_out_edges(rt.graph, node, [&rt, &g, &visited_sinks, &visited_nodes] (const RouteTreeEdge &e) -> void {*/
		const auto &child = get_target(rt.graph, branch);
		check_route_tree_internal(rt, child, g, visited_sinks, visited_nodes);
	}
}

void check_route_tree(const route_tree_t &rt, const net_t &net, GlobalRRGraph &g)
{
	RouteTreeNode rt_root = route_tree_get_rt_node(rt, net.source.rr_node);
	vector<int> sinks;
	for (const auto &s : net.sinks) {
		sinks.push_back(s.rr_node);
	}
	vector<int> visited_sinks;

	vector<int> visited_nodes;
	check_route_tree_internal(rt, rt_root, g, visited_sinks, visited_nodes);
	int num_rt_nodes = 0;
	for (const auto &rt_node : route_tree_get_nodes(rt)) {
		++num_rt_nodes;
	}
	assert(num_rt_nodes == num_edges(rt.graph)+1);
	vector<int> duplicated_nodes;
	sort(begin(visited_nodes), end(visited_nodes));
	int current = visited_nodes[0];
	bool added_duplicated = false;
	for (int i = 1; i < visited_nodes.size(); ++i) {
		/*zlog_info(delta_log, "visited nodes %d\n", visited_nodes[i]);*/
		if (visited_nodes[i] == current) {
			if (!added_duplicated) {
				duplicated_nodes.push_back(visited_nodes[i]);
				added_duplicated = true;
			}
		} else {
			current = visited_nodes[i];
			added_duplicated = false;
		}
	}
	char buffer[256];
	if (!duplicated_nodes.empty()) {
		zlog_error(delta_log, "Error: net %d visited_nodes has duplicates: \n", net.vpr_id);
		for (const auto &d : duplicated_nodes) {
			sprintf_rr_node(d, buffer);
			zlog_error(delta_log, "%s\n", buffer);
		}
		//write_graph(rt.graph, "duplicate.dot",
				//[&duplicated_nodes, &rt] (RouteTreeNode rt_node) -> string {
				//char s_rr_node[256];
				//char buffer[256];
				//const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);
				//if (find(begin(duplicated_nodes), end(duplicated_nodes), rt_node_p.rr_node) !=
						//end(duplicated_nodes)) {
				//sprintf_rr_node(rt_node_p.rr_node, s_rr_node);
				//sprintf(buffer, "label=\"%s\" style=filled fillcolor=red", s_rr_node);
				//} else {
				//sprintf_rr_node(rt_node_p.rr_node, s_rr_node);
				//sprintf(buffer, "label=\"%s\"", s_rr_node);
				//}
				//return buffer;
				//},
				//[] (const RouteTreeEdge &rt_edge) -> string {
				//return "";
				//},
				//[] (const RouteTreeNode &rt_node) -> bool {
				//return false;
				//});
		assert(false);
	}

	sort(visited_sinks.begin(), visited_sinks.end());
	sort(sinks.begin(), sinks.end());

	assert(!sinks.empty());
	
	if (visited_sinks != sinks) {
		zlog_error(delta_log, "Error: Visited %lu sinks out of %lu sinks of net %d\n", visited_sinks.size(), sinks.size(), net.vpr_id);
		vector<int> only_in_required;
		vector<int> only_in_visited;
		vector<int> sym_difference;
		set_difference(sinks.begin(), sinks.end(), visited_sinks.begin(), visited_sinks.end(), back_inserter(only_in_required));
		set_difference(visited_sinks.begin(), visited_sinks.end(), sinks.begin(), sinks.end(), back_inserter(only_in_visited));
		/*set_symmetric_difference(sinks.begin(), sinks.end(), visited_sinks.begin(), visited_sinks.end(), back_inserter(sym_difference));*/
		/*assert(difference == sym_difference);*/

		//write_graph(rt.graph, "error.dot",
				//[&rt] (const RouteTreeNode &rt_node) -> string {
					//char buffer[256];
					//sprintf_rr_node(get_vertex_props(rt.graph, rt_node).rr_node, buffer);
					//return buffer;
				//},
				//[] (const RouteTreeEdge &rt_edge) -> string {
					//return "";
				//},
				//[&rt] (const RouteTreeNode &rt_node) -> bool {
					//return !get_vertex_props(rt.graph, rt_node).valid;
				//});


		zlog_error(delta_log, "Only in required: ");
		for (const int d : only_in_required) {
			zlog_error(delta_log, "%d ", d);
		}
		zlog_error(delta_log, "\n");

		zlog_error(delta_log, "Only in visited: ");
		for (const int d : only_in_visited) {
			zlog_error(delta_log, "%d ", d);
		}
		zlog_error(delta_log, "\n");

		zlog_error(delta_log, "Visited sinks: ");
		for (const int d : visited_sinks) {
			zlog_error(delta_log, "%d ", d);
		}
		zlog_error(delta_log, "\n");

		assert(false);
	}
}

boolean paralar(struct s_router_opts opts,
		float **net_delay, t_slack * slacks, t_ivec ** clb_opins_used_locally, boolean timing_analysis_enabled)
{
	if (setvbuf(stdout, NULL, _IONBF, 0)) {
		exit(-1);
	}

	init_logging();
	
	GlobalRRGraph g;
	map<std::tuple<int, int, int, int, int>, vector<int>> lookup;
	build_global_graph(g, lookup);

	vector<net_t> nets, global_nets;
	init_nets(nets, global_nets, lookup, opts.bb_factor, false);

	init_sprintf_rr_node(&g);

	vector<route_tree_t> route_trees(nets.size());
	for (auto &rt : route_trees) {
		route_tree_init(rt, &g);
	}

	for (const auto &ei : get_edges(g)) {
		auto &edge = get_edge_props(g, ei);
		edge.mult = 0;
	}

	SinkRouter sink_router(g);

	for (int iter = 0; iter < opts.max_router_iterations; ++iter) {
		for (auto &rt : route_trees) {
			route_tree_clear(rt);
		}

		for (const auto &ei : get_edges(g)) {
			auto &edge = get_edge_props(g, ei);
			edge.util = 0;
		}

		for (const auto &net : nets) {
			for (const auto &sink : net.sinks) {
				extern int nx, ny;
				vector<box> bbs { bg::make<box>(0, 0, nx+2, ny+2) };
				float delay;
				vector<GlobalRRNode> added;
				sink_router.route(&net.source, &sink, bbs, 0, route_trees[net.local_id], added, delay, nullptr);
			}
		}

		for (const auto &net : nets) {
			check_route_tree(route_trees[net.local_id], net, g);
		}

		for (const auto &rt : route_trees) {
			for (const auto &rt_edge_i : get_edges(rt.graph)) {
				const auto &rt_edge = get_edge_props(rt.graph, rt_edge_i);
				auto &rr_edge = get_edge_props(g, rt_edge.rr_edge);
				++rr_edge.util;
			}
		}

		int max_util = 0;
		for (const auto &ei : get_edges(g)) {
			auto &edge = get_edge_props(g, ei);
			max_util = std::max(max_util, edge.util);
		}

		const float rho = 0.0001;	
		float norm_T1T2 = 0;
		for (const auto &ei : get_edges(g)) {
			const auto &edge = get_edge_props(g, ei);

			int capacity = edge.ignore ? std::numeric_limits<int>::max() : opts.fixed_channel_width/4;

			int m;
			if (edge.util > capacity) {
				m = 1;
			} else {
				m = 0;
			}
			float gradient = std::max(0.0f, (float) (edge.util - capacity));
			float delay = 1;
			float T1 = delay + m*(edge.mult + rho*gradient);
			float T2 = -gradient;
			norm_T1T2 += sqrt(T1*T1 + T2*T2);	
		}
		float step_size = (1.01/(iter+1))/(norm_T1T2);

		for (const auto &ei : get_edges(g)) {
			auto &edge = get_edge_props(g, ei);

			int capacity = edge.ignore ? std::numeric_limits<int>::max() : opts.fixed_channel_width/4;

			edge.mult += step_size*(std::max(0.0f, (float) (edge.util-capacity)));;
		}

		printf("Iter %d Max util %d Step size %g\n", iter, max_util, step_size);
	}

	return FALSE;
}
