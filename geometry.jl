using Combinatorics, Graphs
using StatsBase, Random, Distributions
using DelimitedFiles, GraphIO
using Memoization

@memoize function splitEdges(n)#See findPercolation.jl for comments on first two functions
    grd=Graphs.grid([n,n],periodic=true)
    gre=unique(map(x->sort(x),collect(Iterators.flatten(map(x->collect(powerset(x, 2, 2)),filter(x->length(x)==4,simplecycles_limited_length(grd,4,4*n^2)))))))
    return [filter(x->has_edge(grd,x[1],x[2]),gre),filter(x->!has_edge(grd,x[1],x[2]),gre)]
end

function ensambleMean(func,n,k,β,μ,default)
    out=Threads.Atomic{Float64}(0)
    edges=splitEdges(n)
    (p1,p2)=(1-1/(1+exp(β*(μ))),1-1/(1+exp(β*(μ-1))))
    Threads.@threads for i=1:k
        (n1,n2)=(rand(Binomial(length(edges[1]),p1)),rand(Binomial(length(edges[2]),p2)))
        if n1+n2>0
            tmp=Graphs.SimpleGraphs.SimpleGraphFromIterator(map(x->Graphs.SimpleEdge(x[1],x[2]),vcat(StatsBase.sample(edges[1],n1,replace=false),StatsBase.sample(edges[2],n2,replace=false))))
            add_vertices!(tmp, n^2-nv(tmp))
            Threads.atomic_add!(out,func(tmp))
        else
            Threads.atomic_add!(out,default)
        end
    end
    return out[]/(k)
end
const n=50
const k=1024
const intervals=50.

t1=@elapsed dt1=[ensambleMean(g->mean(degree(g)),n,k,β,μ,0.) for β=-20:40/intervals:20, μ=-3:6/intervals:3]
writedlm("meanDegree.csv",  dt1, ',')
println("Mean degree took "*string(t1)*"s")

t2=@elapsed dt2=[ensambleMean(g->mean(local_clustering_coefficient(g,vertices(g))),n,k,β,μ,0.) for β=-20:40/intervals:20, μ=-3:6/intervals:3]
writedlm("meanLocalClusteringCoef.csv",  dt2, ',')
println("Local Clustering Coeficcients took "*string(t2)*"s")

println("Total computation time: "*string((t1+t2)/60)*" minutes")