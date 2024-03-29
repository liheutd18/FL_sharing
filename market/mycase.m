function mpc = mycase
%CASE33BW  Power flow data for 33 bus distribution system from Baran & Wu
%    Please see CASEFORMAT for details on the case file format.
%
%    Data from ...
%       M. E. Baran and F. F. Wu, "Network reconfiguration in distribution
%       systems for loss reduction and load balancing," in IEEE Transactions
%       on Power Delivery, vol. 4, no. 2, pp. 1401-1407, Apr 1989.
%       doi: 10.1109/61.25627
%       URL: http://doi.org/10.1109/61.25627

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [  %% (Pd and Qd are specified in kW & kVAr here, converted to MW & MVAr below)
	1	3	0	0	0	0	1	1	0	10	1	1	1;
	2	1	-270	-130	0	0	1	1	0	10	1	1.1	0.9;
	3	1	-150	-70	0	0	1	1	0	10	1	1.1	0.9;
	4	1	520	250	0	0	1	1	0	10	1	1.1	0.9;
	5	1	860	420	0	0	1	1	0	10	1	1.1	0.9;
	6	1	670	320	0	0	1	1	0	10	1	1.1	0.9;
	7	1	2490	1200	0	0	1	1	0	10	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	0	0	10	-10	1	10	1	10	-10	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [  %% (r and x specified in ohms here, converted to p.u. below)
	1	2	0.108	0.108	0	0	0	0	0	0	1	-360	360;
	1	3	0.18	0.18	0	0	0	0	0	0	1	-360	360;
	1	4	0.153	0.153	0	0	0	0	0	0	1	-360	360;
	1	5	0.135	0.135	0	0	0	0	0	0	1	-360	360;
	1	6	0.162	0.162	0	0	0	0	0	0	1	-360	360;
	1	7	0.135	0.135	0	0	0	0	0	0	1	-360	360;
];

%% convert branch impedances from Ohms to p.u.
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
Vbase = mpc.bus(1, BASE_KV) * 1e3;      %% in Volts
Sbase = mpc.baseMVA * 1e6;              %% in VA
mpc.branch(:, [BR_R BR_X]) = mpc.branch(:, [BR_R BR_X]) / (Vbase^2 / Sbase);

%% convert loads from kW to MW
mpc.bus(:, [PD, QD]) = mpc.bus(:, [PD, QD]) / 1e3;