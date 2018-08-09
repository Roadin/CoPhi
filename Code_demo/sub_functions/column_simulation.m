function col = column_simulation(d1, group)

distribution_species = 1:1:3;
dist_index = randsample(distribution_species ,1,true,[0.4, 0.3, 0.3]);
% 1 normal distribution
% 2 possion
% 3 gamma

col = zeros( d1, 1);

numb_obs = floor( d1/group);
if dist_index == 1 % normal distribution
   
    mu = normrnd( 0, 1, group);
    sigma =  unifrnd( 1, 100, group);
    
    for i = 1:( group - 1)
        for j = ( ( i-1 )*numb_obs + 1): i*numb_obs
            col( j ) = normrnd( mu( i), sigma( i ));
        end
    end
    
    for j = ((group - 1)*numb_obs + 1 ): d1
        col( j ) = normrnd( mu( group), sigma( group));
    end
    
elseif dist_index == 2 %poisson distribution
  
    %lambda =randsample(1:1:10, group);
    lambda = unifrnd( 1, 10, group);
    
     for i = 1:( group - 1)
        for j = ( ( i-1 )*numb_obs + 1): i*numb_obs
            col( j ) = poissrnd( lambda( i ));
        end
    end
    
    for j = ((group - 1)*numb_obs + 1 ): d1
        col( j ) =  poissrnd( group );
    end
    
else %gamma distribution
   
     mu =  unifrnd( 1, 10, group);
    sigma =  unifrnd( 1, 10, group);
    
    for i = 1:( group - 1)
        for j = ( ( i-1 )*numb_obs + 1): i*numb_obs
            col( j ) = gamrnd( mu( i), sigma( i ));
        end
    end
    
    for j = ((group - 1)*numb_obs + 1 ): d1
        col( j ) = gamrnd( mu( group), sigma( group));
    end
end

end
    