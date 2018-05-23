using Roots

include("surveHelper.jl")
include("surveSolve.jl")

# Modify supply side to add tax

function F(p)
    return p + log(1+exp(-p))
end

function Fvec(p)
    return p + log.(1+exp.(-p))
end

# Supply function (given price, cost matrix, and counts of cars)
function supplytax(p::Array{Float64,2}, CYZ, my)
    Sz = zeros(10,10);
    for ix in 1:10
        for jx in 1:10
            cyz = exp.(Fvec(p) - CYZ[ix,jx,:,:])
            s = cyz ./ (1 + sum(cyz))

            Sz += my[ix,jx] * s
        end
    end

    return Sz
end

# Supply of y=(i0,j0) for z=(i1,j1)
function syztax(p, CYZ, i0, j0, i1, j1)
    denom = 1 + sum(exp.(Fvec(p) - CYZ[i0,j0,:,:]))
    num = exp(F(p[i1,j1]) - CYZ[i0,j0,i1,j1])
    return num / denom
end

function Wgtax(p, CYZ, my, UXZ, nx)
    E = supplytax(p, CYZ, my) - demand(p, UXZ, nx)
    return E
end

# Gradient at z=(i,j), pz
function WgCoordtax(i, j, pz, p, CYZ, my, UXZ, nx)
    p0 = copy(p);
    p0[i,j] = pz;

    return Wgtax(p0, CYZ, my, UXZ, nx)[i,j]
end

function coordDescentJacobitax(CYZ,my,UXZ,nx, precision)

    pCurr = ones(10,10)
    pNext = ones(10,10)
    itCD = 0

    tic()
    for it in 1:10000

        pNext = copy(pCurr)

        for i in 1:10
            for j in 1:10
                f = pz -> WgCoordtax(i,j,pz,pCurr,CYZ,my,UXZ,nx)
                pNext[i,j] = fzero(f, 0., 2.)
            end
        end

        pCurr = copy(pNext)

        precCurr = prec(supplytax(pCurr,CYZ,my), demand(pCurr,UXZ,nx))
        if precCurr < precision
            itCD = it
            break
        end
    end
    timeCD = toc()

    pCD = pCurr;

    SCD = supplytax(pCD, CYZ, my)
    DCD = demand(pCD, UXZ, nx)
    precCD = prec(SCD, DCD)

    trCD = totalRides(SCD)
    avgprcCD = averagePrice(SCD,pCD)

    CDstats = [precCD, trCD, avgprcCD, itCD, timeCD]
    return CDstats
end


function coordDescentGaussSeideltax(CYZ,my,UXZ,nx, precision)

    pCurr = ones(10,10)
    pNext = ones(10,10)
    itCD = 0

    tic()
    for it in 1:10000

        pNext = copy(pCurr)

        for i in 1:10
            for j in 1:10
                f = pz -> WgCoordtax(i,j,pz,pNext,CYZ,my,UXZ,nx)
                pNext[i,j] = fzero(f, 0., 2.)
            end
        end

        pCurr = copy(pNext)

        precCurr = prec(supplytax(pCurr,CYZ,my), demand(pCurr,UXZ,nx))
        if precCurr < precision
            itCD = it
            break
        end
    end
    timeCD = toc()

    pCD = pCurr;

    SCD = supplytax(pCD, CYZ, my)
    DCD = demand(pCD, UXZ, nx)
    precCD = prec(SCD, DCD)

    trCD = totalRides(SCD)
    avgprcCD = averagePrice(SCD,pCD)

    CDstats = [precCD, trCD, avgprcCD, itCD, timeCD]
    return CDstats, pCD
end
