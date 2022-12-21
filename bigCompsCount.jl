using Combinatorics, Graphs
using StatsBase, Random, Distributions
using DelimitedFiles, GraphIO
using Memoization

@memoize function splitEdges(n)
    grd=Graphs.grid([n,n])
    gre=unique(map(x->sort(x),collect(Iterators.flatten(map(x->collect(powerset(x, 2, 2)),filter(x->length(x)==4,simplecycles_limited_length(grd,4,4*n^2)))))))
    return [filter(x->has_edge(grd,x[1],x[2]),gre),filter(x->!has_edge(grd,x[1],x[2]),gre)]
end

function ensambleMean(func,n,k,β,μ)
    out=Threads.Atomic{Float64}(0)
    edges=splitEdges(n)
    (p1,p2)=(1-1/(1+exp(β*(μ))),1-1/(1+exp(β*(μ-1))))
    Threads.@threads for i=1:k
        (n1,n2)=(rand(Binomial(length(edges[1]),p1)),rand(Binomial(length(edges[2]),p2)))
        if n1+n2>0
            Threads.atomic_add!(out,func(Graphs.SimpleGraphs.SimpleGraphFromIterator(map(x->Graphs.SimpleEdge(x[1],x[2]),vcat(StatsBase.sample(edges[1],n1,replace=false),StatsBase.sample(edges[2],n2,replace=false))))))
        end
    end
    return out[]/(k)
end

function numBigComps(g,cutoff)
    return float(length(filter(x->x>cutoff,map(x->length(x),connected_components(g)))))
end

const n1=50
const n2=100
const n3=200
const k=32
const intervals=50.

t1=@elapsed dt1=[ensambleMean(g->numBigComps(g,n1^2/10),n1,k,β,μ) for β=-20:40/intervals:20, μ=-3:6/intervals:3]
writedlm( "numBigComps50.csv",  dt1, ',')
println("Big Components 50 took "*string(t1)*"s")

t2=@elapsed dt2=[ensambleMean(g->numBigComps(g,n2^2/10),n2,k,β,μ) for β=-20:40/intervals:20, μ=-3:6/intervals:3]
writedlm( "numBigComps100.csv",  dt2, ',')
println("Big Components 100 took "*string(t2)*"s")

t3=@elapsed dt3=[ensambleMean(g->numBigComps(g,n3^2/10),n3,k,β,μ) for β=-20:40/intervals:20, μ=-3:6/intervals:3]
writedlm( "numBigComps200.csv",  dt3, ',')
println("Big Components 200 took "*string(t3)*"s")


println("Total computation time: "*string((t1+t2+t3)/60)*" minutes")
