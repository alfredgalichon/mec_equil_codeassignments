># Aggregate Deferred Acceptance
# (c) 2018 Erica Moszkowski

using Distances, CSV

# For a given type of rider, figure out which type of driver to propose to
# using greedy algorithm

# μxyAx is the set of available matches available to type x
# αxyx is a 2-dimensional matrix containing the preferences of type x
# nxx is a scalar
function propose_one(μxyAx, αxyx, nxx)
    μxyPx              = zeros(10,10);      # Which ys does x propose to?
    total_proposals    = 0.;
    max_iter           = 100;
    t                  = 0;

    while t <= max_iter && total_proposals < nxx

        # Find best counterparty available to x and the utility they generate
        αmax, ymax = findmax(αxyx);
        y1, y2     = ind2sub(size(αxyx), ymax);

        # Take the most of them that you can
        μxyPx[y1,y2] = min(μxyAx[y1,y2], nxx - total_proposals);

        # Update cumulative number of proposals
        total_proposals += μxyPx[y1,y2];

        # Set the utility from that type to -Inf so you don't propose to them again
        αxyx[y1,y2] = -Inf;

        # Iterate
        t += 1;
    end

    return μxyPx;
end

# Figure out proposals for 1 step of algorithm
function propose_all(μxyA,αxy,nx)
    μxyP = zeros(10,10,10,10);
    for x1 in 1:10
        for x2 in 1:10
            μxyP[x1,x2,:,:] = propose_one(μxyA[x1,x2,:,:], αxy[x1,x2,:,:], nx[x1,x2]);
        end
    end
    μxyP;
end

# For a single driver, figure out which proposals to accept
function engagement_one(μxyPy,γxyy,myy)

    μxyEx              = zeros(10,10);      # Which xs does y accept?
    total_accepts      = 0.;
    max_iter           = 100;
    t                  = 0;

    while t <= max_iter && total_accepts < myy

        # Find best counterparty available to y and the utility they generate
        γmax, xmax = findmax(γxyx);
        x1, x2     = ind2sub(size(γxyx), xmax);

        # Take the most of them that you can
        μxyEx[x1,x2] = min(μxyPx[x1,x2], myy - total_accepts);

        # Update cumulative number of proposals
        total_accepts += μxyEx[x1,x2];

        # Set the utility from that type to -Inf so you don't propose to them again
        αxyy[x1,x2] = -Inf;

        # Iterate
        t += 1;
    end

    return μxyEx;
end

# Figure out acceptances for 1 step of algorithm
function engagement_all(μxyP,γxy,my)
    μxyE = zeros(10,10,10,10);
    for y1 in 1:10
        for y2 in 1:10
            μxyE[:,:,y1,y2] = propose_one(μxyP[:,:,y1,y2], γxy[:,:,y1,y2], my[y1,y2]);
        end
    end
    μxyE
end

# Run the algorithm
function aggregate_deferred_acceptance(μA, αxy, γxy, nx, my; max_iter = 1000)

    μAt  = μA;
    μxyP = zeros(10,10,10,10);
    μxyE = zeros(10,10,10,10);

    i = 0
    tic()
    for i = 1:1000
        # Compute proposals
        μxyP = propose_all(μAt, αxy, nx);

        # Compute acceptances
        μxyE = engagement_all(μxyP,γxy,my);

        # Update available offers
        μAt = μAt - (μxyP - μxyE);

        if maximum(abs.(μxyP - μxyE)) < 1e-2;
            break
        end
    end
    run_time = toc()

    return μxyE, i, run_time;
end

# Welfare of riders of type x
function rider_welfare_one(αxyx, μxyx, nxx)
    utility_matched    = sum(μxyx .* αxyx)
    n_matched          = sum(μxyx)
    n_unmatched        = nxx - n_matched
    utility_unmatched  = -2. * n_unmatched

    utility_matched + utility_unmatched
end

# Welfare of drivers of type y
function driver_welfare_one(γxyy, μxyy, myy)
    utility_matched    = sum(μxyy .* γxyy)
    n_matched          = sum(μxyy)
    n_unmatched        = myy - n_matched
    utility_unmatched  = -2. * n_unmatched

    utility_matched + utility_unmatched
end

# Aggregate welfare of riders
function rider_welfare_all(αxy, μxy, nx)
    rider_welfare = zeros(10,10)
    for x1 in 1:10
        for x2 in 1:10
            rider_welfare[x1,x2] = rider_welfare_one(αxy[x1,x2,:,:], μxy[x1,x2,:,:], nx[x1,x2]);
        end
    end
    sum(rider_welfare)
end

# Aggregate welfare of drivers
function rider_welfare_all(γxy, μxy, my)
    driver_welfare = zeros(10,10)
    for y1 in 1:10
        for y2 in 1:10
            driver_welfare[y1,y2] = driver_welfare_one(γxy[:,:,y1,y2], μxy[:,:,y1,y2], my[y1,y2]);
        end
    end
    sum(driver_welfare)
end



##########################################################3
## Setup
##########################################################3

# αxy(x,y) = αxy([x1, x2], [y1,y2]) = αxy(x1, x2, y1, y2)
αxy = zeros(10,10,10,10); # Utility of riders
γxy = zeros(10,10,10,10); # Utility of drivers
nx = zeros(10,10)       ; # Distribution of riders
my = zeros(10,10)       ; # Distribution of drivers
μA = ones(10,10,10,10)  ; # Available matches
μxy = zeros(10,10,10,10) ; # Assigned matches

# Initialize
for x1 = 1:10
    for x2 = 1:10

        # Distribution of riders
        nx[x1,x2] = 1.;

        # Distribution of drivers
        my[x1, x2] = x1/10 * x2/10;

        for y1 = 1:10
            for y2 = 1:10
                # Utility of riders
                αxy[x1,x2,y2,y1] = if euclidean([x1,x2]./10, [y1,y2]./10) >= 0.5
                    -1.;
                else
                    0.;
                end

                # Utility of drivers
                γxy[x1,x2,y1,y2] = -euclidean([x1,x2]./10, [y1,y2]./10);
            end
        end
    end
end

##########################################################
# Run
##########################################################
matches, n_iter, run_time = aggregate_deferred_acceptance(μA, αxy,
                                                          γxy, nx, my);
rider_welfare             = rider_welfare_all(αxy, matches, nx);
driver_welfare            = rider_welfare_all(γxy, matches, my);

results = DataFrame()
results[:rider_welfare]  = rider_welfare
results[:driver_welfare] = driver_welfare
results[:n_iterations]   = n_iter
results[:runtime]        = run_time
CSV.write("aggregate_da_results.csv", results)
