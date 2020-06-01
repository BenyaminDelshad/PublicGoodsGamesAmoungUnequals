% This function calculate the payoff based on threshold func. currently but
% you can change it based on your need. the goal is to provide essential
% details for equilibrium analysis.

function [payoff]= PayoffCalculator(e, x, r, t)
% e is the endowment dist. and x is players contribution rate vector.
% r is the productivity vector and t is the threshold.
% payoff is the n-size vector which is payoff of each players.
    sz = size(e);
    nPlayer = sz(2);
    payoff = zeros(1,nPlayer);
    yvec = e .* r; % players' absolute contributions
    cur_payoff = 0;
    overall_coop_rate = sum(x);
    overal_coop = dot(e,x);
    
    overall_coop_rate = overal_coop; % there are two methods to cosider 
    % as a overall cooperation rate. If you want to use total rate of 
    % cooperation instead of total endowments wich is contirubted, do 
    % comment this line.
    
    %threshold_T = 1/3;
    %threshold_T = 1/2;
    %threshold_T = 2/3;
    threshold_T = t;
    delta = 10; % fixed parameter you can change it. power of sigmoid!
    our_way_of_payoff = 3; % choose your method here!
    % we have three way of having payoff for players.
      
    if (our_way_of_payoff == 1) % way 1: Normally
        cur_payoff = x*yvec';
    
    elseif (our_way_of_payoff == 2) % Way 2: Exact threshold
        if(overall_coop_rate >= threshold_T)
            cur_payoff = x*yvec';
        else
            cur_payoff = 0;
        end
        
    elseif (our_way_of_payoff == 3) % Way 3: Sigmoid functions instead of exact threshold
        cur_payoff = (x*yvec') / (1 + exp(-1 * delta * (overall_coop_rate - threshold_T)));
    end
    
    % Assign payoffs: our payoff divided equally between players.
    payoff(1)=cur_payoff/nPlayer+(1-x(1))*e(1);
    payoff(2)=cur_payoff/nPlayer+(1-x(2))*e(2); 
    payoff(3)=cur_payoff/nPlayer+(1-x(3))*e(3);

end