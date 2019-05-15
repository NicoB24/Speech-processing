function [logProb_alpha,alpha_matrix,logProb_beta,beta_matrix,gamma_matrix,xi_matrix] = logfwd2(x,means,vars,transitions)

% LOGFWD Log version of the forward procedure
%
%    LOGPROB = LOGFWD(X,MEANS,VARS,TRANSITIONS) returns the likelihood of
%    the 2-dimensional sequence X (one observation per row) with respect to
%    a Markov model with N states having means MEANS and variances VARS
%    (stored in N elements lists with empty matrices as first and last
%    elements to symbolize the entry and exit states) and transition matrix
%    TRANSITIONS.
%      Alternately, LOGFWD(X,HMM) can be used, where HMM is an object of the
%    form:
%       HMM.means = MEANS;
%       HMM.vars = VARS;
%       HMM.trans = TRANSITIONS
%
TOLERANCIA = 1e-6; % Tolerancia de error

if nargin == 2,
  model = means;
  means = model.means;
  vars = model.vars;
  model.trans(model.trans<1e-100) = 1e-100;
  logTrans = log(model.trans);
end;

transitions(transitions<1e-100)= 1e-100;
logTrans = log(transitions);
numStates = length(means);
nMinOne = numStates - 1;
[numPts,dim] = size(x);


%Alpha
log2pi = log(2*pi);
for i=2:nMinOne,
  invSig{i} = inv(vars{i});
  logDetVars2(i) = - 0.5 * log(det(vars{i})) - log2pi;
end;

% Initialize the alpha vector for the emitting states
for i=2:nMinOne,
  X = x(1,:)-means{i}';
  alpha(i) = logTrans(1,i) ...
      - 0.5 * (X * invSig{i}) * X' + logDetVars2(i);
end;
alpha = alpha(:);

% Do the forward recursion
alpha_matrix=[];
alpha_matrix=[alpha_matrix alpha];
for t = 2:numPts,
  alphaBefore = alpha;
  for i = 2:nMinOne,
    X = x(t,:)-means{i}';
    alpha(i) = logsum( alphaBefore(2:nMinOne) + logTrans(2:nMinOne,i) ) ...
	- 0.5 * (X * invSig{i}) * X' + logDetVars2(i);
  end;
  alpha_matrix=[alpha_matrix alpha(:)];
end;

%Beta
% Initialize the beta vector for the emitting states
for i=2:nMinOne,
  beta(i) = logTrans(i,numStates);
end;
beta = beta(:);

log2pi = log(2*pi);
for i=2:nMinOne,
  invSig{i} = inv(vars{i});
  logDetVars2(i) = - 0.5 * log(det(vars{i})) - log2pi;
end;

% Do the backward recursion
beta_matrix=[];
beta_matrix=[beta beta_matrix];
for t = numPts:-1:2,
  betaAfter = beta;
  for j=2:nMinOne
    X = x(t,:)-means{j}';
    b(j)=- 0.5 * (X * invSig{j}) * X' + logDetVars2(j);
  end
  for i = 2:nMinOne
    beta(i) = logsum( betaAfter(2:nMinOne) + transpose(logTrans(i,2:nMinOne))+ transpose(b(2:nMinOne)));
  end;
  beta_matrix=[beta(:) beta_matrix];
end;

b= zeros(1,nMinOne);
for j=2:nMinOne
    X = x(1,:)-means{j}';
    b(j)=- 0.5 * (X * invSig{j}) * X' + logDetVars2(j);
end

% Terminate the recursion with the final state
logProb_alpha =  logsum( alpha(2:nMinOne) + logTrans(2:nMinOne,numStates) );
logProb_beta = logsum( b(2:nMinOne) + beta(2:nMinOne,1)' + logTrans(1,2:nMinOne));

%Calculate gamma_matrix
[N,T] = size(alpha_matrix);
gamma_matrix=zeros(N,T);
for i=1:N
    for j=1:T
        gamma_matrix(i,j)=alpha_matrix(i,j)+beta_matrix(i,j)-logsum(alpha_matrix(2:N,i)+beta_matrix(2:N,i));
    end
end

%Calculate xi_matrix
log2pi = log(2*pi);
for i=2:nMinOne,
  invSig{i} = inv(vars{i});
  logDetVars2(i) = - 0.5 * log(det(vars{i})) - log2pi;
end;

alpha=alpha_matrix(2:end,:);
beta=beta_matrix(2:end,:);

logTrans = logTrans(2:end-1,2:end-1);
for t = 2:size(alpha,2)
    divisor = logsum(alpha(:,t) + beta(:,t));
    for k = 1:size(alpha,1)
        X = x(t,:)-means{k+1}';
        b = - 0.5 * (X * invSig{k+1}) * X' + logDetVars2(k+1);        
        for j = 1:size(alpha,1)
            xi_matrix(j,k,t) = alpha(j,t-1) + beta(k,t) + logTrans(j,k) + b - divisor;
        end
    end
end


%=================================
function result = logsum(logv)

len = length(logv);
if (len<2);
  error('Subroutine logsum cannot sum less than 2 terms.');
end;

% First two terms
if (logv(2)<logv(1)),
  result = logv(1) + log( 1 + exp( logv(2)-logv(1) ) );
else,
  result = logv(2) + log( 1 + exp( logv(1)-logv(2) ) );
end;

% Remaining terms
for (i=3:len),
  term = logv(i);
  if (result<term),
    result = term   + log( 1 + exp( result-term ) );
  else,
    result = result + log( 1 + exp( term-result ) );
  end;    
end;
