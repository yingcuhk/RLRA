function mpc = case2(gamma)
%BASED ON																				
%CASE9    Power flow data for 9 bus, 3 generator case.																					
%   Please see CASEFORMAT for details on the case file format.																					
%																					
%   Based on data from Joe H. Chow's book, p. 70.																					
																					
%   MATPOWER																					
%   $Id: case9.m,v 1.11 2010/03/10 18:08:14 ray Exp $																					
																					
%% MATPOWER Case Format : Version 2																					
mpc.version = '2';																					
																					
%%-----  Power Flow Data  -----%%																					
%% system MVA base																					
mpc.baseMVA = 100;																					
																					
%% bus data																					
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin								
mpc.bus = [																					
	1	2	75*gamma	-84.7*gamma	0	0	1	1	0	345	1	1.1	0.9;								
	2	2	105*gamma	22.8*gamma	0	0	1	1	0	345	1	1.1	0.9;								
];																					
																					
%% generator data																					
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [																					
	1	0	0	300	-30	1	100	1	250	75	0	0	0	0	0	0	0	0	0	0	0;
	2	163	0	300	-30	1	100	1	300	70	0	0	0	0	0	0	0	0	0	0	0;
];																					
					
impendence = 0.01008 + 0.0504i;
%% branch data																					
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax								
mpc.branch = [																					
	1	2	real(impendence)	imag(impendence)	0	9900	250	250	0	0	1	-360	360;									
];																					
																					
%%-----  OPF Data  -----%%																					
%% area data																					
%	area	refbus																			
mpc.areas = [																					
	1	2;																			
];																					
																					
%% generator cost data																					
%	1	startup	shutdown	n	x1	y1	...	xn	yn												
%	2	startup	shutdown	n	c(n-1)	...	c0														
mpc.gencost = [																					
	2	1500	0	3	0.11	5	150;														
	2	2000	0	3	0.085	1.2	600;																											
];																					
																					
																					
