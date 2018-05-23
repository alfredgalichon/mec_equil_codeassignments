using Distributions

###############
# Model Setup #
###############

# Set up grids

# Passengers
nx = ones(10,10);

# Drivers
my = zeros(10,10);
for i in 1:10
    for j in 1:10
        my[i,j] = (i/10) * (j/10)
    end
end

# Utility of passengers
function α(x::Array{Int,1},z::Array{Int,1})
    return -1*(norm((x-z)/10,2) >= 0.5)
end

# Compute Preference Terms

# α(x,y) = α(xi,xj,yi,yj)
αxy = zeros(10,10,10,10);

# γ(x,y) = γ(xi,xj,yi,yj)
γxy = zeros(10,10,10,10);

for i0 in 1:10
    for j0 in 1:10
        for i1 in 1:10
            for j1 in 1:10
                αxy[i0,j0,i1,j1] = α([i0,j0],[i1,j1])
                γxy[i0,j0,i1,j1] = -norm([i0,j0]-[i1,j1],2) / 10.
            end
        end
    end
end

################
# Gale-Shapley #
################

# For each x, find μxy (10x10) such that
# 1. μxyP <= μxyA
# 2. sum(μxyP) <= nx

xi = 1
xj = 1

function proposalObjective(μxyP, μxyA, αxy, nx)

    # Check proposal is valid
    if sum(μxyP .<= μxyA) < 100
        return 1e6
    end

    # Check μ is positive
    if sum(μxyP .< 0) > 0
        return 1e6
    end
    
    # Check feasible count
    nMatched = sum(μxyP)
    if nMatched > nx
        return 1e6
    end

    # Compute utility
    matchU = sum(μxyP .* αxy)
    unmatchU = -2*(nx - sum(μxyP))
    return -(matchU + unmatchU)
end

function imax(a)
    i,j = ind2sub(size(a), indmax(a));
    return (i,j)
end

# proposal for a single x (xi, xj)
function proposalGreedy(μxyAx, αxyx, nxx)
    μxyPx = zeros(10,10);
    cumQty = 0.;
    nIt = 0; # cut off if proposed to everyone

    tried = zeros(10,10);
    
    while (cumQty < nxx) & (nIt < 100)
        
        # Identify best counterparty (not already considered)
        iB, jB = imax(αxyx - 10*tried);
        
        # Take maximal amount
        μxyPx[iB,jB] = min(μxyAx[iB,jB], nxx-cumQty);

        # Update cumulative
        cumQty += μxyPx[iB,jB]

        # Update tried
        tried[iB,jB] = 1.;
        
        nIt += 1;
    end

    return μxyPx
end

function proposalIteration(μxyA, αxy, nx)
    μxyP = zeros(10,10,10,10)
    
    for xi in 1:10
        for xj in 1:10
            μxyP[xi,xj,:,:] = proposalGreedy(μxyA[xi,xj,:,:], αxy[xi,xj,:,:], nx[xi,xj])
        end
    end
    
    return μxyP
end


# engagement for a single y (yi, yj)
function engagementGreedy(μxyPy, γxyy, myy)
    μxyEy = zeros(10,10);
    cumQty = 0.;
    nIt = 0; # cut off if proposed to everyone

    tried = zeros(10,10);
    
    while (cumQty < myy) & (nIt < 100)
        # Identify best counterparty
        iB, jB = imax(γxyy - 10*tried);
        
        # Take maximal amount
        μxyEy[iB,jB] = min(μxyPy[iB,jB], myy-cumQty);

        # Update cumulative
        cumQty += μxyEy[iB,jB]

        # Update tried
        tried[iB,jB] = 1.;
        
        nIt += 1;
    end

    return μxyEy
end

function engagementIteration(μxyP, γxy, my)
    μxyE = zeros(10,10,10,10)
    
    for yi in 1:10
        for yj in 1:10
            μxyE[:,:,yi,yj] = engagementGreedy(μxyP[:,:,yi,yj], γxy[:,:,yi,yj], my[yi,yj])
        end
    end
    
    return μxyE
end

# Run Gale-Shapley

μxyA0 = ones(10,10,10,10);
μxyA1 = ones(10,10,10,10);
μxyP1 = ones(10,10,10,10);
μxyE1 = ones(10,10,10,10);

GSit = 0;

tic()
for it in 1:1000

    μxyP1 = proposalIteration(μxyA0, αxy, nx);
    μxyE1 = engagementIteration(μxyP1, γxy, my);
    μxyA1 = μxyA1 - (μxyP1 - μxyE1);

    # Check converge
    if sum(abs.(μxyP1 - μxyE1)) < 1e-2
        print("Converged")
        GSit = it;
        break
    end
    
    μxyA0 = copy(μxyA1);
end
GStime = toc()

@show GStime, GSit

# Compute summary stats

function passWelfare(μxy, αxy, nx)
    uMatch = sum(μxy .* αxy);
    uNomatch = -2 * sum(nx - sum(μxy,[3,4])[:,:,1,1]);

    return uMatch + uNomatch
end

function carWelfare(μxy, γxy, my)
    uMatch = sum(μxy .* γxy);
    uNomatch = -2 * sum(my - sum(μxy,[1,2])[1,1,:,:]);

    return uMatch + uNomatch
end

@show GSpass = passWelfare(μxyP1, αxy, nx)
@show GScar = carWelfare(μxyP1, γxy, my)

######################
# Adachi's Algorithm #
######################


# General tiebreaker noise for α and γ
srand(1);
αε = rand(Uniform(), 10, 10, 10, 10)*0.01;
γε = rand(Uniform(), 10, 10, 10, 10)*0.01;

αxyε = αxy + αε;
γxyε = γxy + γε;


function updateU(αxyε, γxyε, vj)
    ui = zeros(10,10);
    for xi in 1:10
        for yi in 1:10
            # γij >= vj 
            eligibleY = γxyε[xi,yi,:,:] .>= vj
            ui[xi,yi] = max(maximum(αxyε[xi,yi,:,:] -10*(1 - eligibleY)), -2.)
        end
    end
    return ui
end

function updateV(αxyε, γxyε, ui)
    vj = zeros(10,10);
    for xj in 1:10
        for yj in 1:10
            # αij >= ui 
            eligibleX = αxyε[:,:,xj,yj] .>= ui
            vj[xj,yj] = max(maximum(γxyε[:,:,xj,yj] -10*(1 - eligibleX)), -2.)
        end
    end
    return vj
end


# Inital point
ui0 = maximum(αxyε,[3,4])[:,:,1,1];
vj0 = ones(10,10)*-2;
ui1 = copy(ui0)
vj1 = copy(vj0)

Ait = 0

tic()
for it in 1:100
    ui1 = updateU(αxyε, γxyε, vj0)
    vj1 = updateV(αxyε, γxyε, ui1)

    @show it
    @show prec = norm(ui1 - ui0) + norm(vj1 - vj0)
    @show norm(ui1-ui0)
    if prec < 1e-5
        print("CONVERGED")
        Ait = it
        break
    end

    ui0 = copy(ui1)
    vj0 = copy(vj1)
end
Atime = toc()

@show Apass = sum(ui1)
@show Acar = sum(vj1)

#######################
# GS with tiebreakers #
#######################

# Run Gale-Shapley

μxyA0 = ones(10,10,10,10);
μxyA1 = ones(10,10,10,10);
μxyP1 = ones(10,10,10,10);
μxyE1 = ones(10,10,10,10);

GSεit = 0;

tic()
for it in 1:1000

    μxyP1 = proposalIteration(μxyA0, αxyε, nx);
    μxyE1 = engagementIteration(μxyP1, γxyε, nx);
    μxyA1 = μxyA1 - (μxyP1 - μxyE1);

    # Check converge
    if sum(abs.(μxyP1 - μxyE1)) < 1e-2
        print("Converged")
        GSεit = it;
        break
    end
    
    μxyA0 = copy(μxyA1);
end
GSεtime = toc()

uiGS = sum(μxyP1 .* αxyε, [3,4])[:,:,1,1]
vjGS = sum(μxyP1 .* γxyε, [1,2])[1,1,:,:]

@show GSεpass = passWelfare(μxyP1, αxyε, nx)
@show GSεcar = carWelfare(μxyP1, γxyε, nx)

################
# Save to file #
################

using DataFrames
using CSV

df = DataFrame(hcat([GSit, GStime, GSpass, GScar], [Ait, Atime, Apass, Acar],
                    [GSεit, GSεtime, GSεpass, GSεcar]))
CSV.write("table.csv", df)
