using Roots

# Code to optimize using each method and helpers

include("surgeHelper.jl")

####################
# gradient descent #
####################

# Gradient of objective W = excess supply function
# so ∇W = S - D
function Wg(p, CYZ, my, UXZ, nx)
    E = supply(p, CYZ, my) - demand(p, UXZ, nx)
    return E
end

function gradientDescent(CYZ, my, UXZ, nx, precision)

    pCurr = ones(10,10)
    pNext = ones(10,10)
    ε = 0.2
    itSGD = 0
    
    tic()
    for it in 1:10000
        pNext = pCurr - ε*Wg(pCurr, CYZ, my, UXZ, nx)
        pCurr = pNext

        precCurr = prec(supply(pCurr,CYZ,my), demand(pCurr,UXZ,nx))
        if precCurr < precision
            itSGD = it
            break
        end
    end
    timeSGD = toc()

    pGD = pCurr;

    SGD = supply(pGD, CYZ, my)
    DGD = demand(pGD, UXZ, nx)
    precSGD = prec(SGD, DGD)

    trSGD = totalRides(SGD)
    avgprcSGD = averagePrice(SGD,pGD)

    GDstats = [precSGD, trSGD, avgprcSGD, itSGD, timeSGD]
    return GDstats
end

##################
# Newton Descent #
##################


# Hessian of objective W
function WH(p, CYZ, my, UXZ, nx)
    H = zeros(100,100);

    # Compute μyz and μxz terms
    SYZ = zeros(10,10,10,10);
    DXZ = zeros(10,10,10,10);
    for i0 in 1:10
        for j0 in 1:10
            SYZ[i0,j0,:,:] = reshape([syz(p, CYZ, i0, j0, iz, jz) for iz in 1:10 for jz in 1:10],10,10);
            DXZ[i0,j0,:,:] = reshape([dxz(p, UXZ, i0, j0, iz, jz) for iz in 1:10 for jz in 1:10],10,10);
        end
    end
    
    # derivative of W(i1,j1) wrt (i0,j0)
    for i0 in 1:10
        for j0 in 1:10

            # supply by each y for z=(i0,j0)
            s0 = SYZ[:,:,i0,j0]
            d0 = DXZ[:,:,i0,j0]

            for i1 in 1:10
                for j1 in 1:10

                    s1 = SYZ[:,:,i1,j1]
                    d1 = DXZ[:,:,i1,j1]
                    
                    if (i0 == i1) & (j0==j1)
                        H[(i0-1)*10+j0,(i0-1)*10+j0] = sum(my .* s0 .* (1-s0)) - sum(nx .* d0 .* (1-d0))
                    else
                        H[(i0-1)*10+j0,(i1-1)*10+j1] = sum(my .* s0 .* s1) - sum(nx .* d0 .* d1)
                    end
                end
            end
        end
    end

    #H = H + (H - diagm(diag(H)))';
    return H
end

function newtonMethod(CYZ, my, UXZ, nx, precision)

    pCurr = ones(100,1)*2
    pNext = copy(pCurr)
    ε = 0.2
    itND = 0

    tic()
    for it in 1:1000

        Hess = WH(reshape(pCurr,10,10), CYZ, my, UXZ, nx);
        Grad = reshape(Wg(reshape(pCurr,10,10), CYZ, my, UXZ, nx),100,1);
        increment = -ε*inv(Hess)*Grad;
        pNext = pCurr + increment;
        pCurr = pNext;
        #reshape(pCurr,10,10)

        precCurr = prec(supply(reshape(pCurr,10,10),CYZ,my), demand(reshape(pCurr,10,10),UXZ,nx))
        if precCurr < precision
            itND = it
            break
        end
    end
    timeND = toc()

    pND = reshape(pCurr,10,10);

    SND = supply(pND, CYZ, my)
    DND = demand(pND, UXZ, nx)
    precND = prec(SND, DND)

    trND = totalRides(SND)
    avgprcND = averagePrice(SND,pND)

    NDstats = [precND, trND, avgprcND, itND, timeND]
    return NDstats
end


###############################
# Coordinate Descent (Jacobi) #
###############################


# Gradient at z=(i,j), pz
function WgCoord(i, j, pz, p, CYZ, my, UXZ, nx)
    p0 = copy(p);
    p0[i,j] = pz;

    return Wg(p0, CYZ, my, UXZ, nx)[i,j]
end

function coordDescentJacobi(CYZ,my,UXZ,nx,precision)

    pCurr = ones(10,10)
    pNext = ones(10,10)
    itCD = 0

    tic()
    for it in 1:10000

        pNext = copy(pCurr)

        for i in 1:10
            for j in 1:10
                f = pz -> WgCoord(i,j,pz,pCurr,CYZ,my,UXZ,nx)
                pNext[i,j] = fzero(f, 0., 2.)
            end
        end

        pCurr = copy(pNext)

        precCurr = prec(supply(pCurr,CYZ,my), demand(pCurr,UXZ,nx))
        if precCurr < precision
            itCD = it
            break
        end
    end
    timeCD = toc()

    pCD = pCurr;

    SCD = supply(pCD, CYZ, my)
    DCD = demand(pCD, UXZ, nx)
    precCD = prec(SCD, DCD)

    trCD = totalRides(SCD)
    avgprcCD = averagePrice(SCD,pCD)

    CDstats = [precCD, trCD, avgprcCD, itCD, timeCD]
    return CDstats
end



#####################################
# Coordinate Descent (Gauss-Seidel) #
#####################################


function coordDescentGaussSeidel(CYZ,my,UXZ,nx,precision)

    pCurr = ones(10,10)
    pNext = ones(10,10)
    itCD = 0

    tic()
    for it in 1:10000

        pNext = copy(pCurr)

        for i in 1:10
            for j in 1:10
                f = pz -> WgCoord(i,j,pz,pNext,CYZ,my,UXZ,nx)
                pNext[i,j] = fzero(f, 0., 2.)
            end
        end

        pCurr = copy(pNext)

        precCurr = prec(supply(pCurr,CYZ,my), demand(pCurr,UXZ,nx))
        if precCurr < precision
            itCD = it
            break
        end
    end
    timeCD = toc()

    pCD = pCurr;

    SCD = supply(pCD, CYZ, my)
    DCD = demand(pCD, UXZ, nx)
    precCD = prec(SCD, DCD)

    trCD = totalRides(SCD)
    avgprcCD = averagePrice(SCD,pCD)

    CDstats = [precCD, trCD, avgprcCD, itCD, timeCD]
    return CDstats
end
