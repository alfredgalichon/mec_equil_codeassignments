# Adachi's algorithm
# (c) 2018 Erica Moszkowski

using Distances, Distributions, CSV

##########################################################3
## Functions
##########################################################3

# αxyx is a 2x2 matrix of the utility x gets from each y
# γxyx is a 2x2 matrix of the utility each y gets from x
# nxx is the number of type-x agents there are
# vy is a 2x2 matrix of the utility that y is currently guaranteed
function update_rider(αxyx, γxyx, nxx, vy)

    # Figure out who x can offer a better deal to
    available_y_inds = γxyx .>= vy

    # Make a matrix of the utilities rider x gets from
    # the ys available to them (all other ys give -Inf)
    ux_available = fill(-Inf, 10, 10);
    ux_available[available_y_inds] = αxyx[available_y_inds];

    # Find the y that gives the maximum utility
    umax, xmatch = findmax(ux_available)
    umax_all = max(umax, -2.)
    return umax_all
end

# αxyy is a 2x2 matrix of the utility each x gets from y
# γxyy is a 2x2 matrix of the utility y gets from each x
# myy is the number of type-x agents there are
# ux is a 2x2 matrix of the utility that x is currently guaranteed
function update_driver(γxyy, αxyy, myy, ux)

    # Figure out who x can offer a better deal to
    available_x_inds = αxyy .>= ux

    # Make a matrix of the utilities driver y gets from
    # the xs available to them (all other xs give utility -Inf)
    vy_available = fill(-Inf, 10, 10);
    vy_available[available_x_inds] = γxyy[available_x_inds];

    # Find the x that gives the maximum utility
    vmax, ymatch = findmax(vy_available)
    vmax_all     = max(vmax, -2.)
    return vmax_all
end

function update_riders(αxy, γxy, nx, vy)
    ux                = zeros(10,10);
    for x1 in 1:10
        for x2 in 1:10
            ux[x1, x2] = update_rider(αxy[x1,x2,:,:], γxy[x1,x2,:,:], nx[x1,x2], vy)
        end
    end
    ux
end

function update_drivers(γxy, αxy, my, ux)
    vy            = zeros(10,10);
    for y1 in 1:10
        for y2 in 1:10
            vy[y1,y2] = update_driver(γxy[:,:,y1,y2], αxy[:,:,y1,y2], my[y1,y2], ux)
        end
    end
    vy
end


function adachi(αxy, γxy, nx, my)
    # Initialize ux at the utility their best match would give to x
    ux_old = zeros(10,10);
    for x1 in 1:10
        for x2 in 1:10
            ux_old[x1,x2] = maximum(αxy[x1,x2,:,:]);
        end
    end
    ux_current = ux_old;

    # Initialize py at outside option utility
    vy_old     = fill(-2., 10, 10)
    vy_current = vy_old;

    # Iterate
    i = 0
    tic()
    for i = 1:1000
        # Update
        ux_current = update_riders(αxy,γxy,nx,vy_old)
        vy_current = update_drivers(γxy,αxy,my,ux_current)

        # Check if we're done
        uxdiff = abs.(ux_current - ux_old);
        vydiff = abs.(vy_current - vy_old);
        if maximum(uxdiff) .<= 1e-2 && maximum(vydiff) .<= 1e-2
            break
        end

        # Increment
        ux_old = ux_current
        vy_old = vy_current
    end
    run_time = toc()

    return ux_current, vy_current, i, run_time
end


##########################################################3
## Setup
##########################################################3

# αxy(x,y) = αxy([x1, x2], [y1,y2]) = αxy(x1, x2, y1, y2)
αxy = zeros(10,10,10,10); # Utility of riders
γxy = zeros(10,10,10,10); # Utility of drivers
nx = zeros(10,10)       ; # Distribution of riders
my = zeros(10,10)       ; # Distribution of drivers

# Initialize
for x1 = 1:10
    for x2 = 1:10

        # Distribution of riders
        nx[x1,x2] = 1.;

        # Distribution of drivers
        # my[x1, x2] = x1/10 * x2/10;
        my[x1, x2] = 1.;

        for y1 = 1:10
            for y2 = 1:10
                # Utility of riders
                αxy[x1,x2,y1,y2] = if euclidean([x1,x2]./10, [y1,y2]./10) >= 0.5
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

# Add random noise
srand(0)
unif = Uniform(0., 0.1)
αxy = αxy + rand(unif, size(αxy)...);
γxy = γxy + rand(unif, size(γxy)...);

##########################################################
## Run
##########################################################
ux, vy, n_iter, run_time = adachi(αxy, γxy, nx, my);
rider_welfare = sum(ux)
driver_welfare = sum(vy)

results = DataFrame()
results[:rider_welfare]  = rider_welfare
results[:driver_welfare] = driver_welfare
results[:n_iterations]   = n_iterations
results[:runtime]        = run_time

CSV.write("adachi_results.csv", results)

# jldopen("adachi_out.jld", "w") do file
#     write(file, "results", results)
# end
