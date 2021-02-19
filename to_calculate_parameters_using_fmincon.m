function [ fitting_results,size_population_raw,x,x2,xS,x3,simple_,double_,simple_stretched_,triple_ ] = to_calculate_parameters_using_fmincon( binranges,data,models )

% fitting using fmincon function

%   X = FMINCON(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS) minimizes with
%   the default optimization parameters replaced by values in OPTIONS, an
%   argument created with the OPTIMOPTIONS function. See OPTIMOPTIONS for
%   details. For a list of options accepted by FMINCON refer to the
%   documentation.

%   X = FMINCON(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure with the function FUN in PROBLEM.objective, the start point
%   in PROBLEM.x0, the linear inequality constraints in PROBLEM.Aineq
%   and PROBLEM.bineq, the linear equality constraints in PROBLEM.Aeq and
%   PROBLEM.beq, the lower bounds in PROBLEM.lb, the upper bounds in
%   PROBLEM.ub, the nonlinear constraint function in PROBLEM.nonlcon, the
%   options structure in PROBLEM.options, and solver name 'fmincon' in
%   PROBLEM.solver. Use this syntax to solve at the command line a problem
%   exported from OPTIMTOOL. The structure PROBLEM must have all the fields.
%
%   [X,FVAL] = FMINCON(FUN,X0,...) returns the value of the objective
%   function FUN at the solution X.
%
%   [X,FVAL,EXITFLAG] = FMINCON(FUN,X0,...) returns an EXITFLAG that
%   describes the exit condition of FMINCON. Possible values of EXITFLAG
%   and the corresponding exit conditions are listed below. See the
%   documentation for a complete description.
%
%   All algorithms:
%     1  First order optimality conditions satisfied.
%     0  Too many function evaluations or iterations.
%    -1  Stopped by output/plot function.
%    -2  No feasible point found.
%   Trust-region-reflective, interior-point, and sqp:
%     2  Change in X too small.
%   Trust-region-reflective:
%     3  Change in objective function too small.
%   Active-set only:
%     4  Computed search direction too small.
%     5  Predicted change in objective function too small.
%   Interior-point and sqp:
%    -3  Problem seems unbounded.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = FMINCON(FUN,X0,...) returns a structure
%   OUTPUT with information such as total number of iterations, and final
%   objective function value. See the documentation for a complete list.
%
%   [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = FMINCON(FUN,X0,...) returns the
%   Lagrange multipliers at the solution X: LAMBDA.lower for LB,
%   LAMBDA.upper for UB, LAMBDA.ineqlin is for the linear inequalities,
%   LAMBDA.eqlin is for the linear equalities, LAMBDA.ineqnonlin is for the
%   nonlinear inequalities, and LAMBDA.eqnonlin is for the nonlinear
%   equalities.
%
%   [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD] = FMINCON(FUN,X0,...) returns the
%   value of the gradient of FUN at the solution X.
%
%   [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = FMINCON(FUN,X0,...)
%   returns the value of the exact or approximate Hessian of the Lagrangian
%   at X.

size_population_raw = nansum(data(:,2));
fitting_results.size_population = size_population_raw;


%Options=optimset('TolFun',1e-6,'TolX',1e-15,'MaxFunEvals',1000, 'MaxIter' , 1000);
%Options=optimset('TolFun',1e-6,'TolX',1e-15,'MaxFunEvals',1000, 'MaxIter' , 1000, 'GradObj', 'on');
%Options=optimset('TolFun',1e-6,'TolX',1e-15,'MaxFunEvals',1000, 'MaxIter' , 1000, 'GradObj', 'on', 'algorithm', 'trust-region-reflective', 'Hessian', 'user-supplied');
%Options=optimset('TolFun',1e-6,'TolX',1e-15,'MaxFunEvals',1000, 'MaxIter' , 1000, 'GradObj', 'on', 'algorithm', 'interior-point', 'Hessian', 'user-supplied','HessFcn','@hessinterior');
%Options=optimset('TolFun',1e-15,'TolX',1e-15);
Options=optimset('TolFun',1e-6,'TolX',1e-15,'MaxFunEvals',1000, 'MaxIter' , 1000, 'GradObj', 'on', 'GradConstr', 'on');

fitting_results.options = Options;


if ~isempty(find(ismember(models,'MonoExpo'),1))
    
    %[x, ~, exitflag]=fmincon(@(par) Loglikelihood2_deriv(par, data,'MonoExpo', @simple_exp2_beta), [size_population_raw/100 1],...
    %    [], [], [], [],[0 0], [size_population_raw/2 10],[],Options);
    
    % try to add non linear constraints
    [x, ~, exitflag]=fmincon(@(par) Loglikelihood2_deriv(par, data,'MonoExpo', @simple_exp2_beta), [size_population_raw/100 1],...
       [], [], [], [],[0 0.1], [size_population_raw/2 10],@(par) calc_pop_simple(par,data,size_population_raw),Options);
    
    fitting_results.MonoExpo.parameters = x;
    fitting_results.MonoExpo.nparam = 2;
    fitting_results.MonoExpo.A = x(1);
    fitting_results.MonoExpo.T = 1/x(2);
    fitting_results.MonoExpo.flag = exitflag;
    size_population = floor(vpa(sum(simple_exp2_beta(x,data)),3));
    fitting_results.MonoExpo.size_population = double(size_population);
    simple_=simple_exp2_beta(x,data);
    residues = (simple_(:) - data(:,2)) ./data(:,2);
    residues(~isfinite(residues))=NaN;
    fitting_results.MonoExpo.residuals = cat(1,binranges,transpose(residues));
    clear exitflag
    clear size_population
    clear residues
else
    x = [];
    simple_ = [];
end

if ~isempty(find(ismember(models,'DoubleExpo'),1))
    
    %[x2, ~, exitflag]=fmincon(@(par) Loglikelihood2_deriv(par, data,'DoubleExpo', @double_exp2_beta), [size_population_raw/100 0.5 size_population_raw/100 2], ...
    %    [], [], [], [],[0 0 0 0], [size_population_raw/2 10 size_population_raw/2 10],[],Options);
    
    % try to add non linear constraints
    [x2, ~, exitflag]=fmincon(@(par) Loglikelihood2_deriv(par, data,'DoubleExpo', @double_exp2_beta), [size_population_raw/100 0.5 size_population_raw/100 2], ...
        [], [], [], [],[0 0.1 0 0.1], [size_population_raw/2 10 size_population_raw/2 10],@(par) calc_pop_double(par,data,size_population_raw),Options);
    
    fitting_results.DoubleExpo.parameters = x2;
    fitting_results.DoubleExpo.nparam = 4;
    if 1/x2(2) < 1/x2(4)
        fitting_results.DoubleExpo.B1 = x2(1);
        fitting_results.DoubleExpo.T1 = 1/x2(2);
        fitting_results.DoubleExpo.B2 = x2(3);
        fitting_results.DoubleExpo.T2 = 1/x2(4);
        fitting_results.DoubleExpo.flag = exitflag;
        size_population = floor(vpa(sum(double_exp2_beta(x2,data)),3));
        fitting_results.DoubleExpo.size_population=double(size_population);
        size_population1 = floor(vpa(sum(simple_exp2_beta([x2(1) x2(2)],data)),3));
        fitting_results.DoubleExpo.size_population1=double(size_population1);
        fitting_results.DoubleExpo.percent_population1=double(size_population1)/double(size_population)*100;
        size_population2 = floor(vpa(sum(simple_exp2_beta([x2(3) x2(4)],data)),3));
        fitting_results.DoubleExpo.size_population2=double(size_population2);
        fitting_results.DoubleExpo.percent_population2=double(size_population2)/double(size_population)*100;
        double_=double_exp2_beta(x2,data);
        residues = (double_(:) - data(:,2)) ./data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.DoubleExpo.residuals = cat(1,binranges,transpose(residues));
    elseif 1/x2(4) < 1/x2(2)
        fitting_results.DoubleExpo.B1 = x2(3);
        fitting_results.DoubleExpo.T1 = 1/x2(4);
        fitting_results.DoubleExpo.B2 = x2(1);
        fitting_results.DoubleExpo.T2 = 1/x2(2);
        fitting_results.DoubleExpo.flag = exitflag;
        size_population = floor(vpa(sum(double_exp2_beta(x2,data)),3));
        fitting_results.DoubleExpo.size_population=double(size_population);
        size_population1 = floor(vpa(sum(simple_exp2_beta([x2(3) x2(4)],data)),3));
        fitting_results.DoubleExpo.size_population1=double(size_population1);
        fitting_results.DoubleExpo.percent_population1=double(size_population1)/double(size_population)*100;
        size_population2 = floor(vpa(sum(simple_exp2_beta([x2(1) x2(2)],data)),3));
        fitting_results.DoubleExpo.size_population2=double(size_population2);
        fitting_results.DoubleExpo.percent_population2=double(size_population2)/double(size_population)*100;
        double_=double_exp2_beta(x2,data);
        residues = (double_(:) - data(:,2)) ./data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.DoubleExpo.residuals = cat(1,binranges,transpose(residues));
    end
    clear exitflag
    clear size_population size_population1 size_population2
    clear residues
else
    x2 = [];
    double_ = [];
end

if ~isempty(find(ismember(models,'MonoExpo_stretched'),1))
    
    %[xS, ~, exitflag]=fmincon(@(par) Loglikelihood2_deriv(par, data,'MonoExpo_stretched', @simple_exp2_stretched_beta), [size_population_raw/100 1 1],...
    %    [], [], [], [],[0 0 0], [size_population_raw/2 10 10],[],Options);
    
    % try to add non linear constraints
    [xS, ~, exitflag]=fmincon(@(par) Loglikelihood2_deriv(par, data,'MonoExpo_stretched', @simple_exp2_stretched_beta), [size_population_raw/100 1 1],...
        [], [], [], [],[0 0.1 0], [size_population_raw/2 10 10],@(par) calc_pop_simple_stretched(par,data,size_population_raw),Options);
    
    fitting_results.MonoExpo_stretched.parameters = xS;
    fitting_results.MonoExpo_stretched.nparam = 3;
    fitting_results.MonoExpo_stretched.D = xS(1);
    fitting_results.MonoExpo_stretched.T = 1/xS(2);
    fitting_results.MonoExpo_stretched.power = xS(3);
    fitting_results.MonoExpo_stretched.flag = exitflag;
    size_population = floor(vpa(sum(simple_exp2_stretched_beta(xS,data)),3));
    fitting_results.MonoExpo_stretched.size_population=double(size_population);
    simple_stretched_=simple_exp2_stretched_beta(xS,data);
    residues = (simple_stretched_(:) - data(:,2)) ./data(:,2);
    residues(~isfinite(residues))=NaN;
    fitting_results.MonoExpo_stretched.residuals = cat(1,binranges,transpose(residues));
    clear exitflag
    clear size_population
    clear residues
else
    xS = [];
    simple_stretched_ = [];
end

if ~isempty(find(ismember(models,'TripleExpo'),1))
    
    %[x3, ~, exitflag]=fmincon(@(par) Loglikelihood2_deriv(par, data,'TripleExpo', @triple_exp2_beta), [size_population_raw/100 0.5 size_population_raw/100 2 size_population_raw/100 5 ], ...
    %    [], [], [], [],[0 0 0 0 0 0], [size_population_raw/2 10 size_population_raw/2 10 size_population_raw/2 10],[],Options);
    
    % try to add non linear constraints
    [x3, ~, exitflag]=fmincon(@(par) Loglikelihood2_deriv(par, data,'TripleExpo', @triple_exp2_beta), [size_population_raw/100 0.5 size_population_raw/100 2 size_population_raw/100 5 ], ...
        [], [], [], [],[0 0.1 0 0.1 0 0.1], [size_population_raw/2 10 size_population_raw/2 10 size_population_raw/2 10],@(par) calc_pop_triple(par,data,size_population_raw),Options);
    
    fitting_results.TripleExpo.parameters = x3;
    fitting_results.TripleExpo.nparam = 6;
    if 1/x3(2) < 1/x3(4) < 1/x3(6)
        fitting_results.TripleExpo.C1 = x3(1);
        fitting_results.TripleExpo.TT1 = 1/x3(2);
        fitting_results.TripleExpo.C2 = x3(3);
        fitting_results.TripleExpo.TT2 = 1/x3(4);
        fitting_results.TripleExpo.C3 = x3(5);
        fitting_results.TripleExpo.TT3 = 1/x3(6);
        fitting_results.TripleExpo.flag = exitflag;
        size_population = floor(vpa(sum(triple_exp2_beta(x3,data)),3));
        fitting_results.TripleExpo.size_population=double(size_population);
        size_population1 = floor(vpa(sum(simple_exp2_beta([x3(1) x3(2)],data)),3));
        fitting_results.TripleExpo.size_population1=double(size_population1);
        fitting_results.TripleExpo.percent_population1=double(size_population1)/double(size_population)*100;
        size_population2 = floor(vpa(sum(simple_exp2_beta([x3(3) x3(4)],data)),3));
        fitting_results.TripleExpo.size_population2=double(size_population2);
        fitting_results.TripleExpo.percent_population2=double(size_population2)/double(size_population)*100;
        size_population3 = floor(vpa(sum(simple_exp2_beta([x3(5) x3(6)],data)),3));
        fitting_results.TripleExpo.size_population3=double(size_population3);
        fitting_results.TripleExpo.percent_population3=double(size_population3)/double(size_population)*100;
        triple_=triple_exp2_beta(x3,data);
        residues = (triple_(:) - data(:,2)) ./data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.TripleExpo.residuals = cat(1,binranges,transpose(residues));
    elseif 1/x3(2)  < 1/x3(6) < 1/x3(4)
        fitting_results.TripleExpo.C1 = x3(1);
        fitting_results.TripleExpo.TT1 = 1/x3(2);
        fitting_results.TripleExpo.C2 = x3(5);
        fitting_results.TripleExpo.TT2 = 1/x3(6);
        fitting_results.TripleExpo.C3 = x3(3);
        fitting_results.TripleExpo.TT3 = 1/x3(4);
        fitting_results.TripleExpo.flag = exitflag;
        size_population = floor(vpa(sum(triple_exp2_beta(x3,data)),3));
        fitting_results.TripleExpo.size_population=double(size_population);
        size_population1 = floor(vpa(sum(simple_exp2_beta([x3(1) x3(2)],data)),3));
        fitting_results.TripleExpo.size_population1=double(size_population1);
        fitting_results.TripleExpo.percent_population1=double(size_population1)/double(size_population)*100;
        size_population2 = floor(vpa(sum(simple_exp2_beta([x3(5) x3(6)],data)),3));
        fitting_results.TripleExpo.size_population2=double(size_population2);
        fitting_results.TripleExpo.percent_population2=double(size_population2)/double(size_population)*100;
        size_population3 = floor(vpa(sum(simple_exp2_beta([x3(3) x3(4)],data)),3));
        fitting_results.TripleExpo.size_population3=double(size_population3);
        fitting_results.TripleExpo.percent_population3=double(size_population3)/double(size_population)*100;
        triple_=triple_exp2_beta(x3,data);
        residues = (triple_(:) - data(:,2)) ./data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.TripleExpo.residuals = cat(1,binranges,transpose(residues));
    elseif 1/x3(4)  < 1/x3(2) < 1/x3(6)
        fitting_results.TripleExpo.C1 = x3(3);
        fitting_results.TripleExpo.TT1 = 1/x3(4);
        fitting_results.TripleExpo.C2 = x3(1);
        fitting_results.TripleExpo.TT2 = 1/x3(2);
        fitting_results.TripleExpo.C3 = x3(5);
        fitting_results.TripleExpo.TT3 = 1/x3(6);
        fitting_results.TripleExpo.flag = exitflag;
        size_population = floor(vpa(sum(triple_exp2_beta(x3,data)),3));
        fitting_results.TripleExpo.size_population=double(size_population);
        size_population1 = floor(vpa(sum(simple_exp2_beta([x3(3) x3(4)],data)),3));
        fitting_results.TripleExpo.size_population1=double(size_population1);
        fitting_results.TripleExpo.percent_population1=double(size_population1)/double(size_population)*100;
        size_population2 = floor(vpa(sum(simple_exp2_beta([x3(1) x3(2)],data)),3));
        fitting_results.TripleExpo.size_population2=double(size_population2);
        fitting_results.TripleExpo.percent_population2=double(size_population2)/double(size_population)*100;
        size_population3 = floor(vpa(sum(simple_exp2_beta([x3(5) x3(6)],data)),3));
        fitting_results.TripleExpo.size_population3=double(size_population3);
        fitting_results.TripleExpo.percent_population3=double(size_population3)/double(size_population)*100;
        triple_=triple_exp2_beta(x3,data);
        residues = (triple_(:) - data(:,2)) ./data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.TripleExpo.residuals = cat(1,binranges,transpose(residues));
    elseif 1/x3(4)  < 1/x3(6) < 1/x3(2)
        fitting_results.TripleExpo.C1 = x3(3);
        fitting_results.TripleExpo.TT1 = 1/x3(4);
        fitting_results.TripleExpo.C2 = x3(5);
        fitting_results.TripleExpo.TT2 = 1/x3(6);
        fitting_results.TripleExpo.C3 = x3(1);
        fitting_results.TripleExpo.TT3 = 1/x3(2);
        fitting_results.TripleExpo.flag = exitflag;
        size_population = floor(vpa(sum(triple_exp2_beta(x3,data)),3));
        fitting_results.TripleExpo.size_population=double(size_population);
        size_population1 = floor(vpa(sum(simple_exp2_beta([x3(3) x3(4)],data)),3));
        fitting_results.TripleExpo.size_population1=double(size_population1);
        fitting_results.TripleExpo.percent_population1=double(size_population1)/double(size_population)*100;
        size_population2 = floor(vpa(sum(simple_exp2_beta([x3(5) x3(6)],data)),3));
        fitting_results.TripleExpo.size_population2=double(size_population2);
        fitting_results.TripleExpo.percent_population2=double(size_population2)/double(size_population)*100;
        size_population3 = floor(vpa(sum(simple_exp2_beta([x3(1) x3(2)],data)),3));
        fitting_results.TripleExpo.size_population3=double(size_population3);
        fitting_results.TripleExpo.percent_population3=double(size_population3)/double(size_population)*100;
        triple_=triple_exp2_beta(x3,data);
        residues = (triple_(:) - data(:,2)) ./data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.TripleExpo.residuals = cat(1,binranges,transpose(residues));
    elseif 1/x3(6)  < 1/x3(2) < 1/x3(4)
        fitting_results.TripleExpo.C1 = x3(5);
        fitting_results.TripleExpo.TT1 = 1/x3(6);
        fitting_results.TripleExpo.C2 = x3(1);
        fitting_results.TripleExpo.TT2 = 1/x3(2);
        fitting_results.TripleExpo.C3 = x3(3);
        fitting_results.TripleExpo.TT3 = 1/x3(4);
        fitting_results.TripleExpo.flag = exitflag;
        size_population = floor(vpa(sum(triple_exp2_beta(x3,data)),3));
        fitting_results.TripleExpo.size_population=double(size_population);
        size_population1 = floor(vpa(sum(simple_exp2_beta([x3(5) x3(6)],data)),3));
        fitting_results.TripleExpo.size_population1=double(size_population1);
        fitting_results.TripleExpo.percent_population1=double(size_population1)/double(size_population)*100;
        size_population2 = floor(vpa(sum(simple_exp2_beta([x3(1) x3(2)],data)),3));
        fitting_results.TripleExpo.size_population2=double(size_population2);
        fitting_results.TripleExpo.percent_population2=double(size_population2)/double(size_population)*100;
        size_population3 = floor(vpa(sum(simple_exp2_beta([x3(3) x3(4)],data)),3));
        fitting_results.TripleExpo.size_population3=double(size_population3);
        fitting_results.TripleExpo.percent_population3=double(size_population3)/double(size_population)*100;
        triple_=triple_exp2_beta(x3,data);
        residues = (triple_(:) - data(:,2)) ./data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.TripleExpo.residuals = cat(1,binranges,transpose(residues));
    elseif 1/x3(6)  < 1/x3(4) < 1/x3(2)
        fitting_results.TripleExpo.C1 = x3(5);
        fitting_results.TripleExpo.TT1 = 1/x3(6);
        fitting_results.TripleExpo.C2 = x3(3);
        fitting_results.TripleExpo.TT2 = 1/x3(4);
        fitting_results.TripleExpo.C3 = x3(1);
        fitting_results.TripleExpo.TT3 = 1/x3(2);
        fitting_results.TripleExpo.flag = exitflag;
        size_population = floor(vpa(sum(triple_exp2_beta(x3,data)),3));
        fitting_results.TripleExpo.size_population=double(size_population);
        size_population1 = floor(vpa(sum(simple_exp2_beta([x3(5) x3(6)],data)),3));
        fitting_results.TripleExpo.size_population1=double(size_population1);
        fitting_results.TripleExpo.percent_population1=double(size_population1)/double(size_population)*100;
        size_population2 = floor(vpa(sum(simple_exp2_beta([x3(3) x3(4)],data)),3));
        fitting_results.TripleExpo.size_population2=double(size_population2);
        fitting_results.TripleExpo.percent_population2=double(size_population2)/double(size_population)*100;
        size_population3 = floor(vpa(sum(simple_exp2_beta([x3(1) x3(2)],data)),3));
        fitting_results.TripleExpo.size_population3=double(size_population3);
        fitting_results.TripleExpo.percent_population3=double(size_population3)/double(size_population)*100;
        triple_=triple_exp2_beta(x3,data);
        residues = (triple_(:) - data(:,2)) ./data(:,2);
        residues(~isfinite(residues))=NaN;
        fitting_results.TripleExpo.residuals = cat(1,binranges,transpose(residues));
    end
    clear exitflag
    clear size_population size_population1 size_population2 size_population3
    clear residues
else
    x3 = [];
    triple_ = [];
end

% notice before residues are ponderated ones!!

end

