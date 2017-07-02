#ifndef PARALAR_H
#define PARALAR_H

boolean paralar(struct s_router_opts router_opts,
		float **net_delay, t_slack * slacks, t_ivec ** clb_opins_used_locally, boolean timing_analysis_enabled);

#endif
