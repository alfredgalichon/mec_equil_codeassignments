using Roots

# Code to set up model and compute summary statistics

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
function u(x::Array{Int,1},z::Array{Int,1})
    return -10*(norm((x-z)/10,2) >= 0.1)
end

# Compute all cross-x-u values
# U(x,z)
UXZ = zeros(10,10,10,10);

for i0 in 1:10
    for j0 in 1:10
        for i1 in 1:10
            for j1 in 1:10
                UXZ[i0,j0,i1,j1] = u([i0,j0],[i1,j1])
            end
        end
    end
end

# Compute all cross-y-z values
# C(y,z)
CYZ = zeros(10,10,10,10);

for i0 in 1:10
    for j0 in 1:10
        for i1 in 1:10
            for j1 in 1:10
                CYZ[i0,j0,i1,j1] = norm([i0,j0]-[i1,j1],2)^2 / 100.
            end
        end
    end
end


# Demand function (given prices, utility matrix, and counts of passengers)
function demand(p::Array{Float64,2}, UXZ, nx)
    Dz = zeros(10,10);
    for ix in 1:10
        for jx in 1:10
            uxz = exp.(UXZ[ix,jx,:,:] - p)
            d = uxz ./ (1 + sum(uxz))

            Dz += nx[ix,jx] * d 
        end
    end

    return Dz
end

# Supply function (given price, cost matrix, and counts of cars)
function supply(p::Array{Float64,2}, CYZ, my)
    Sz = zeros(10,10);
    for ix in 1:10
        for jx in 1:10
            cyz = exp.(p - CYZ[ix,jx,:,:])
            s = cyz ./ (1 + sum(cyz))

            Sz += my[ix,jx] * s
        end
    end

    return Sz
end

# Demand of x=(i0,j0) for z=(i1,j1)
function dxz(p, CYZ, i0, j0, i1, j1)
    denom = 1 + sum(exp.(UXZ[i0,j0,:,:] - p))
    num = exp(UXZ[i0,j0,i1,j1] - p[i1,j1])
    return num / denom
end

# Supply of y=(i0,j0) for z=(i1,j1)
function syz(p, CYZ, i0, j0, i1, j1)
    denom = 1 + sum(exp.(p - CYZ[i0,j0,:,:]))
    num = exp(p[i1,j1] - CYZ[i0,j0,i1,j1])
    return num / denom
end


#############################
# Set of summary statistics #
#############################

function prec(S, D)
    return maximum(abs.(S - D) ./ (S + D))
end

function totalRides(S)
    return sum(S)
end

function averagePrice(S,P)
    return sum(S .* P) / sum(S)
end

