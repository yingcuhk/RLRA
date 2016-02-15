function mpc = get_mpc(om)
%GET_MPC  Returns the MATPOWER case struct.
%   MPC = GET_MPC(OM)
%
%   See also OPT_MODEL.

%   MATPOWER
%   Copyright (c) 2008-2015 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   $Id: get_mpc.m 2644 2015-03-11 19:34:22Z ray $
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

mpc = om.mpc;
