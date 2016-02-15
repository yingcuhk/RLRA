function rv = mpver(varargin)
%MPVER  Prints or returns MATPOWER version info for current installation.
%   V = MPVER returns the current MATPOWER version number.
%   V = MPVER('all') returns a struct with the fields Name, Version,
%   Release and Date (all strings). Calling MPVER without assigning the
%   return value prints the version and release date of the current
%   installation of MATPOWER, MATLAB (or Octave), the Optimization Toolbox,
%   MIPS and any optional MATPOWER packages.

%   MATPOWER
%   Copyright (c) 2005-2015 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   $Id: mpver.m 2664 2015-03-20 21:31:02Z ray $
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% the following only works if MATPOWER is explicitly in the path,
%% but not if it is only in the current working directory
% fs = filesep;
% p = fileparts(which('runpf'));
% if ~strcmp(p(1),fs)
%   [t, p] = strtok(p, filesep);
% end
% p = p(2:end);
% v{1} = ver(p);

v{1} = struct(  'Name',     'MATPOWER', ... 
                'Version',  '5.1', ...
                'Release',  '', ...
                'Date',     '20-Mar-2015' );
if nargout > 0
    if nargin > 0
        rv = v{1};
    else
        rv = v{1}.Version;
    end
else
    if have_fcn('octave')
        v{2} = ver('octave');
    else
        v{2} = ver('matlab');
    end
    v{3} = ver('optim');
    for n = 1:3
        if n == 3 && isempty(v{3})
            fprintf('\n%-22s -- not installed --', 'Optimization Toolbox');
            continue;
        elseif n == 3 && ~license('test', 'optimization_toolbox')
            fprintf('\n%-22s -- no license --', 'Optimization Toolbox');
            continue;
        end
        fprintf('\n%-22s Version %-9s', v{n}.Name, v{n}.Version);
        if ~isempty(v{n}.Date)
            fprintf('  %11s', v{n}.Date);
            if ~isempty(v{n}.Release)
                fprintf('   Release: %-10s', v{n}.Release);
            end
        end
    end
    fprintf('\n');
    mipsver;
    if have_fcn('sdp_pf')
        sdp_pf_ver;
    else
        fprintf('%-22s -- not installed --\n', 'SDP_PF');
    end
    if have_fcn('yalmip')
        s = have_fcn('yalmip', 'all');
        fprintf('%-22s Version %-10s %-11s\n', 'YALMIP', s.vstr, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'YALMIP');
    end
    if have_fcn('bpmpd')
        if exist('bpver', 'file') == 2
            bpver;
        else
            fprintf('%-22s Version 2.21 or earlier\n', 'BPMPD_MEX');
        end
    else
        fprintf('%-22s -- not installed --\n', 'BPMPD_MEX');
    end
    if have_fcn('clp')
        s = have_fcn('clp', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'CLP', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'CLP');
    end
    if have_fcn('cplex')
        s = have_fcn('cplex', 'all');
        fprintf('%-22s Version %-10s %-11s\n', 'CPLEX', s.vstr, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'CPLEX');
    end
    if have_fcn('glpk')
        s = have_fcn('glpk', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'GLPK', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'GLPK');
    end
    if have_fcn('gurobi')
        gurobiver;
    else
        fprintf('%-22s -- not installed --\n', 'Gurobi');
    end
    if have_fcn('ipopt')
        s = have_fcn('ipopt', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'IPOPT', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'IPOPT');
    end
    if have_fcn('knitro')
        s = have_fcn('knitro', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'KNITRO', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'KNITRO');
    end
    if have_fcn('minopf')
        if exist('minopfver', 'file') == 2
            minopfver;
        else
            fprintf('%-22s Version 3.0b2 or earlier\n', 'MINOPF');
        end
    else
        fprintf('%-22s -- not installed --\n', 'MINOPF');
    end
    if have_fcn('mosek')
        s = have_fcn('mosek', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'MOSEK', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'MOSEK');
    end
    if have_fcn('pardiso')
        s = have_fcn('pardiso', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'PARDISO', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'PARDISO');
    end
    if have_fcn('pdipmopf')
        pdipmopfver;
    else
        fprintf('%-22s -- not installed --\n', 'PDIPMOPF');
    end
    if have_fcn('scpdipmopf')
        scpdipmopfver;
    else
        fprintf('%-22s -- not installed --\n', 'SCPDIPMOPF');
    end
    if have_fcn('sdpt3')
        s = have_fcn('sdpt3', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'SDPT3', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'SDPT3');
    end
    if have_fcn('sedumi')
        s = have_fcn('sedumi', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'SeDuMi', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'SeDuMi');
    end
    if have_fcn('tralmopf')
        tralmopfver;
    else
        fprintf('%-22s -- not installed --\n', 'TRALMOPF');
    end

    fprintf('%-22s %s\n\n', 'Architecture:', computer);
    
    fprintf('  MATPOWER %s is distributed under the 3-clause BSD License.\n', v{1}.Version);
    fprintf('  Please see the LICENSE file for details.\n\n');
end
