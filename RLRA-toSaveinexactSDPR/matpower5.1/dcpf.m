function Va = dcpf(B, Pbus, Va0, ref, pv, pq)
%DCPF  Solves a DC power flow.
%   [VA, SUCCESS] = DCPF(B, PBUS, VA0, REF, PV, PQ) solves for the bus
%   voltage angles at all but the reference bus, given the full system
%   B matrix and the vector of bus real power injections, the initial
%   vector of bus voltage angles (in radians), and column vectors with
%   the lists of bus indices for the swing bus, PV buses, and PQ buses,
%   respectively. Returns a vector of bus voltage angles in radians.
%
%   See also RUNDCPF, RUNPF.

%   MATPOWER
%   Copyright (c) 1996-2015 by Power System Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   and Ray Zimmerman, PSERC Cornell
%
%   $Id: dcpf.m 2644 2015-03-11 19:34:22Z ray $
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% initialize result vector
Va = Va0;

%% update angles for non-reference buses
Va([pv; pq]) = B([pv; pq], [pv; pq]) \ (Pbus([pv; pq]) - B([pv; pq], ref) * Va0(ref));
