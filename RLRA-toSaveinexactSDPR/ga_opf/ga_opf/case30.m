function [baseMVA, bus, gen, branch, areas, gencost] = case30
%CASE30    Power flow data for 30 bus, 6 generator case.
%   Please see 'help caseformat' for details on the case file format.
%
%   Based on data from ...
%     Alsac, O. & Stott, B., "Optimal Load Flow with Steady State Security",
%     IEEE Transactions on Power Apparatus and Systems, Vol. PAS 93, No. 3,
%     1974, pp. 745-751.
%   ... with branch parameters rounded to nearest 0.01, shunt values divided
%   by 100 and shunt on bus 10 moved to bus 5, load at bus 5 zeroed out.
%   Generator locations, costs and limits and bus areas were taken from ...
%     Ferrero, R.W., Shahidehpour, S.M., Ramesh, V.C., "Transaction analysis
%     in deregulated power systems using game theory", IEEE Transactions on
%     Power Systems, Vol. 12, No. 3, Aug 1997, pp. 1340-1347.
%   Generator Q limits were derived from Alsac & Stott, using their Pmax
%   capacities. V limits and line |S| limits taken from Alsac & Stott.

%   MATPOWER
%   $Id: case30.m,v 1.7 2004/09/21 02:41:49 ray Exp $

%%-----  Power Flow Data  -----%%
%% system MVA base
baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
bus = [
	1	3	0	0	0	0	1	1	0	135	1	1.05	0.95;
	2	2	21.7	12.7	0	0	1	1	0	135	1	1.1	0.95;
	3	1	2.4	1.2	0	0	1	1	0	135	1	1.05	0.95;
	4	1	7.6	1.6	0	0	1	1	0	135	1	1.05	0.95;
	5	1	0	0	0	0.19	1	1	0	135	1	1.05	0.95;
	6	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	7	1	22.8	10.9	0	0	1	1	0	135	1	1.05	0.95;
	8	1	30	30	0	0	1	1	0	135	1	1.05	0.95;
	9	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	10	1	5.8	2	0	0	3	1	0	135	1	1.05	0.95;
	11	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	12	1	11.2	7.5	0	0	2	1	0	135	1	1.05	0.95;
	13	2	0	0	0	0	2	1	0	135	1	1.1	0.95;
	14	1	6.2	1.6	0	0	2	1	0	135	1	1.05	0.95;
	15	1	8.2	2.5	0	0	2	1	0	135	1	1.05	0.95;
	16	1	3.5	1.8	0	0	2	1	0	135	1	1.05	0.95;
	17	1	9	5.8	0	0	2	1	0	135	1	1.05	0.95;
	18	1	3.2	0.9	0	0	2	1	0	135	1	1.05	0.95;
	19	1	9.5	3.4	0	0	2	1	0	135	1	1.05	0.95;
	20	1	2.2	0.7	0	0	2	1	0	135	1	1.05	0.95;
	21	1	17.5	11.2	0	0	3	1	0	135	1	1.05	0.95;
	22	2	0	0	0	0	3	1	0	135	1	1.1	0.95;
	23	2	3.2	1.6	0	0	2	1	0	135	1	1.1	0.95;
	24	1	8.7	6.7	0	0.04	3	1	0	135	1	1.05	0.95;
	25	1	0	0	0	0	3	1	0	135	1	1.05	0.95;
	26	1	3.5	2.3	0	0	3	1	0	135	1	1.05	0.95;
	27	2	0	0	0	0	3	1	0	135	1	1.1	0.95;
	28	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	29	1	2.4	0.9	0	0	3	1	0	135	1	1.05	0.95;
	30	1	10.6	1.9	0	0	3	1	0	135	1	1.05	0.95;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin
gen = [
	1	23.54	0	150	-20	1	100	1	80	0;
	2	60.97	0	60	-20	1	100	1	80	0;
	22	21.59	0	62.5	-15	1	100	1	50	0;
	27	26.91	0	48.7	-15	1	100	1	55	0;
	23	19.2	0	40	-10	1	100	1	30	0;
	13	37	0	44.7	-15	1	100	1	40	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status
linedata = [
	1	2	0.02	0.06	0.03	130	130	130	0	0	1;
	1	3	0.05	0.19	0.02	130	130	130	0	0	1;
	2	4	0.06	0.17	0.02	65	65	65	0	0	1;
	3	4	0.01	0.04	0	130	130	130	0	0	1;
	2	5	0.05	0.2	0.02	130	130	130	0	0	1;
	2	6	0.06	0.18	0.02	65	65	65	0	0	1;
	4	6	0.01	0.04	0	90	90	90	0	0	1;
	5	7	0.05	0.12	0.01	70	70	70	0	0	1;
	6	7	0.03	0.08	0.01	130	130	130	0	0	1;
	6	8	0.01	0.04	0	32	32	32	0	0	1;
	6	9	0	0.21	0	65	65	65	0	0	1;
	6	10	0	0.56	0	32	32	32	0	0	1;
	9	11	0	0.21	0	65	65	65	0	0	1;
	9	10	0	0.11	0	65	65	65	0	0	1;
	4	12	0	0.26	0	65	65	65	0	0	1;
	12	13	0	0.14	0	65	65	65	0	0	1;
	12	14	0.12	0.26	0	32	32	32	0	0	1;
	12	15	0.07	0.13	0	32	32	32	0	0	1;
	12	16	0.09	0.2	0	32	32	32	0	0	1;
	14	15	0.22	0.2	0	16	16	16	0	0	1;
	16	17	0.08	0.19	0	16	16	16	0	0	1;
	15	18	0.11	0.22	0	16	16	16	0	0	1;
	18	19	0.06	0.13	0	16	16	16	0	0	1;
	19	20	0.03	0.07	0	32	32	32	0	0	1;
	10	20	0.09	0.21	0	32	32	32	0	0	1;
	10	17	0.03	0.08	0	32	32	32	0	0	1;
	10	21	0.03	0.07	0	32	32	32	0	0	1;
	10	22	0.07	0.15	0	32	32	32	0	0	1;
	21	22	0.01	0.02	0	32	32	32	0	0	1;
	15	23	0.1	0.2	0	16	16	16	0	0	1;
	22	24	0.12	0.18	0	16	16	16	0	0	1;
	23	24	0.13	0.27	0	16	16	16	0	0	1;
	24	25	0.19	0.33	0	16	16	16	0	0	1;
	25	26	0.25	0.38	0	16	16	16	0	0	1;
	25	27	0.11	0.21	0	16	16	16	0	0	1;
	28	27	0	0.4	0	65	65	65	0	0	1;
	27	29	0.22	0.42	0	16	16	16	0	0	1;
	27	30	0.32	0.6	0	16	16	16	0	0	1;
	29	30	0.24	0.45	0	16	16	16	0	0	1;
	8	28	0.06	0.2	0.02	32	32	32	0	0	1;
	6	28	0.02	0.06	0.01	32	32	32	0	0	1;
];

%%-----  OPF Data  -----%%
%% area data
areas = [
	1	8;
	2	23;
	3	26;
];

%% generator cost data
%	1	startup	shutdown	n	x0	y0	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
gencost = [
	2	0	0	3	0.02	2	0;
	2	0	0	3	0.0175	1.75	0;
	2	0	0	3	0.0625	1	0;
	2	0	0	3	0.00834	3.25	0;
	2	0	0	3	0.025	3	0;
	2	0	0	3	0.025	3	0;
];

return;
