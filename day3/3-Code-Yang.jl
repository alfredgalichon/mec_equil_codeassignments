using Gurobi
using JuMP
using CSV
using DataFrames

#############
# Load data #
#############

arcs = CSV.read("NYC_subway/arcs_clean.csv",allowmissing=:none)
nodes = CSV.read("NYC_subway/nodes_clean.csv",allowmissing=:none)

##############
# Data Setup #
##############

source = 452;
sink = 471;

fromNb = arcs[:from_stop_nb]
toNb = arcs[:to_stop_nb]
dis = arcs[:dis_line]

stopName = nodes[:stop_name]

nNodes = max(maximum(fromNb), maximum(toNb))
nArcs = length(fromNb)

# Construct da
da = dis;

# Construct s
s = zeros(nNodes);
s[source] = -1;
s[sink] = 1;

# Construct ∇
A = zeros(nNodes, nArcs);
for n in 1:nArcs
    A[fromNb[n], n] = 1
    A[toNb[n], n] = -1
end
∇=A';

#################################
# Set up JuMP model with Gurobi #
#################################

m = Model(solver=GurobiSolver(Presolve=0))

# Add variables (μ ≥ 0)
@variable(m, μ[1:nArcs] >= 0);

# Add objective (min ∑ₐ dₐ μₐ)
@objective(m, :Min, sum(da[i] * μ[i] for i=1:nArcs))

# Add constraint (∇'μ = s)
@constraint(m, [i=1:nNodes], sum(A[i,:] .* μ) == s[i])

# Solve
status = solve(m)

objval = getobjectivevalue(m)
μSoln = getvalue(μ)

##############
# Print Path #
##############

solnMask = μSoln.==1;
solnFrom = fromNb[solnMask]
solnTo = toNb[solnMask]

path = zeros(Int,length(solnFrom)+1)
path[1] = source
path[end] = sink

currNode = source
for i in 2:(length(solnTo))
    path[i] = solnFrom[find(solnTo.==path[i-1])[1]]
end

lkup = x -> stopName[x];
pathName = [lkup(x) for x in path];

@show pathName

################
# Bellman Ford #
################

# Initialize prices
pz = ones(nNodes)*Inf;
pz[source] = 0

pred = zeros(nNodes)*NaN;
steps = zeros(Int,nNodes);

for i in 1:(nNodes-1)
    for a in 1:nArcs
        if pz[fromNb[a]] + da[a] < pz[toNb[a]]
            pz[toNb[a]] = pz[fromNb[a]] + da[a]
            pred[toNb[a]] = fromNb[a]
            steps[toNb[a]] = steps[fromNb[a]]+1
        end
    end
end

# Construct path
pathBF = zeros(Int,steps[sink]+1)
pathBF[1] = sink
for i in 2:length(pathBF)
    pathBF[i] = pred[pathBF[i-1]]
end
pathBF = reverse(pathBF)

########################
# Approximate Min Cost #
########################

# Define objective
function objective(p)
    reg = 0
    for a in 1:nArcs
        reg -= exp(p[fromNb[a]] - p[toNb[a]] - da[a])
    end
    return s' * p + reg
end

function grad(p)
    g = copy(s)
    for a in 1:nArcs
        dev = exp(p[fromNb[a]] - p[toNb[a]] - da[a])
        g[fromNb[a]] -= dev
        g[toNb[a]] += dev
    end
    return g
end

p = zeros(501)
p[source] = -600.
p[sink] = 500.
@show objective(p)

ε=5.
nIt = 0
err = 0

for it in 1:100000
    ε = norm(grad(p))
    p1 = p + ε*grad(p)

    err = norm(p - p1)
    if it % 1000 == 0
        @show it
        @show objective(p1)
        @show norm(grad(p1))
    end
    
    #@show 
    if err < 1e-5
        @show nIt = it
        break
    end

    p = copy(p1)
end
print("NO")

@show err
@show objective(p)
