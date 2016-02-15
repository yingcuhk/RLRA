function [varargout] = dcopf(varargin)
%DCOPF  Solves a DC optimal power flow.
%   This is a simple wrapper function around OPF that sets the 'model'
%   option to 'DC' before calling OPF.
%   See OPF for the details of input and output arguments.
%
%   See also RUNDCOPF.

%   MATPOWER
%   Copyright (c) 1996-2015 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   $Id: dcopf.m 2644 2015-03-11 19:34:22Z ray $
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[mpc, mpopt] = opf_args(varargin{:});
mpopt = mpoption(mpopt, 'model', 'DC');
[varargout{1:nargout}] = opf(mpc, mpopt);
