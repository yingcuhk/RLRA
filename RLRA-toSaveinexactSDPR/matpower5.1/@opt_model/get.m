function val = get(om, varargin)
%GET  Returns the value of a field.
%   VAL = GET(OM, FIELD1, FIELD2, ...)
%
%   Example:
%       var_order = get(om, 'var', 'order');
%
%   See also OPT_MODEL.

%   MATPOWER
%   Copyright (c) 2008-2015 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   $Id: get.m 2644 2015-03-11 19:34:22Z ray $
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

val = om;
for k = 1:length(varargin)
    if ischar(varargin{k})
        val = val.(varargin{k});
    else
        val = val(varargin{k});
    end
end
