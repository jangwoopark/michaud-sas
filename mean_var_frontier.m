function [weights,means,standard_devs]=mean_var_frontier(mu,sigma,m)

%this program constructs the Markowitz mean-variance frontier for a given
%set of parameters (consisting of a column vector of means, a
%variance-covariance matrix, and the number of steps in the partition of
%the mean dimension.

n=length(mu);
sorted_mean=sort(mu);
max_mean=sorted_mean(n,1);

%We will also probably need to increase the maximum number of iterations
%for 'quadprog' (the quadratic programming solver), and in order to do this
%we also must set a starting point for the solver.  Also, because of the 
%quadratic program being solved, large scale optimization is not possible:

x0quadprog=ones(n,1)./n;
options=optimset('MaxIter',Inf,'LargeScale','off');

%find the global minimum variance portfolio using 'quadprog':
[gmv_weight,gmv_var] = quadprog(sigma,zeros(n,1),zeros(n,n),zeros(n,1),ones(1,n),1,zeros(n,1),ones(n,1),x0quadprog,options);
gmv_mean=mu'*gmv_weight;

%create placeholders
weights=zeros(m,n);
means=zeros(m,1);
standard_devs=zeros(m,1);

%construct the Markowitz mean-variance frontier
stepsize_mean=(max_mean-gmv_mean)/(m-1);
Aeq=[mu';ones(1,n)];
beq=ones(2,1);
for i=1:m;
    means(i,1)=gmv_mean+(stepsize_mean*(i-1));
    beq(1,1)=means(i,1);
    [x,fval] = quadprog(sigma,zeros(n,1),zeros(n,n),zeros(n,1),Aeq,beq,zeros(n,1),ones(n,1),x0quadprog,options);
    weights(i,:)=x';
    standard_devs(i,1)=sqrt(2*fval);
end
