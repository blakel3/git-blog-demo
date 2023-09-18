using Distances

const dt=1e-3
const tmax=2
const frameinterval=.125
const saveinterval=frameinterval

const ld=6
const d0=1
const g0=2
const Y=4e6
const zeta=2e2
const noise=.02/dt
const fc=sqrt(d0)*Y

x=Float64[0]
y=Float64[0]
l=Float64[ld/2]
g=Float64[g0]
theta=Float64[0]

xcell=Any[]
ycell=Any[]
gcell=Any[]
lcell=Any[]
thetacell=Any[]
tcell=Float64[]

t::Float64=0
n::Int64=1

while t<tmax
    global t=t+dt

    ## increase cell lengths and split long cells

    global l=l.+g*dt;
    TO_SPLIT=findall(l.>ld)
    l2=copy(l)
    l2[TO_SPLIT]=l[TO_SPLIT]/2 .-d0/2
    l2=append!(l2,l[TO_SPLIT]/2 .-d0/2)
    
    theta2=append!(theta,theta[TO_SPLIT])

    x2 = copy(x)
    x2[TO_SPLIT] = x[TO_SPLIT] .- .25 * cos.(theta[TO_SPLIT]) .* (l[TO_SPLIT] .+ d0)
    global x = append!(x2 , x[TO_SPLIT] .+ .25 * cos.(theta[TO_SPLIT]) .* (l[TO_SPLIT] .+ d0))

    y2 = copy(y)
    y2[TO_SPLIT] = y[TO_SPLIT] .- .25 * sin.(theta[TO_SPLIT]) .* (l[TO_SPLIT] .+ d0)
    global y = append!(y2 , y[TO_SPLIT] .+ .25 * sin.(theta[TO_SPLIT]) .* (l[TO_SPLIT] .+ d0))

    g2 = copy(g)
    g2[TO_SPLIT] = g0 .+ (rand(length(TO_SPLIT)) .- .5) .* g0
    global g = append!(g2 , g0 .+ (rand(length(TO_SPLIT)) .- .5) .* g0)

    global l = l2
    global theta = theta2

    xpoints = hcat(x .- .5 .* l .* cos.(theta), x, x .+ .5 .* l .* cos.(theta))
    ypoints = hcat(y .- .5 .* l .* sin.(theta), y, y .+ .5 .* l .* sin.(theta))

    D = pairwise(Euclidean(), hcat(x,y))

    CLOSE=zeros(Bool,length(x),length(x))
    for i=1:length(x)
        for j=1:length(x)
            CLOSE[i,j]=l[i]/2+l[j]/2+d0 > D[i,j]
            CLOSE[j,i]=CLOSE[i,j]
        end
    end
    ind=findall(CLOSE)
    for i=1:length(ind)
        
    end



end

