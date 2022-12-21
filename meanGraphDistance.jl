using Combinatorics, Graphs
using StatsBase, Random, Distributions
using DelimitedFiles, GraphIO
using Memoization

@memoize function splitEdges(n)#See findPercolation.jl for comments on first two functions
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

function meanDist(g)#Estimate the mean graph distance (inside components)
    verts=StatsBase.sample(vertices(g),100)
    return mean(map(x->mean(filter(x->x<typemax(Int),dijkstra_shortest_paths(g,x).dists)),verts))
end

const n1=50
const n2=100
const n3=200
const n4=400
const n5=800
const k=128
const intervals=50

t1=@elapsed dt1=[ensambleMean(meanDist,n1,k,β,-β) for β=0.5:0.5/intervals:1]*3/(2*n1)
writedlm( "graphDiams50.csv",  dt1, ',')
println("Graph Diameters 50 took "*string(t1)*"s")

t2=@elapsed dt2=[ensambleMean(meanDist,n2,k,β,-β) for β=0.5:0.5/intervals:1]*3/(2*n2)
writedlm( "graphDiams100.csv",  dt2, ',')
println("Graph Diameters 100 took "*string(t2)*"s")

t3=@elapsed dt3=[ensambleMean(meanDist,n3,k,β,-β) for β=0.5:0.5/intervals:1]*3/(2*n3)
writedlm( "graphDiams200.csv",  dt3, ',')
println("Graph Diameters 200 took "*string(t3)*"s")

t4=@elapsed dt4=[ensambleMean(meanDist,n4,k,β,-β) for β=0.5:0.5/intervals:1]*3/(2*n4)
writedlm( "graphDiams400.csv",  dt4, ',')
println("Graph Diameters 400 took "*string(t4)*"s")

t5=@elapsed dt5=[ensambleMean(meanDist,n5,k,β,-β) for β=0.5:0.5/intervals:1]*3/(2*n5)
writedlm( "graphDiams800.csv",  dt5, ',')
println("Graph Diameters 800 took "*string(t5)*"s")


println("Total computation time: "*string((t1+t2+t3+t4+t5)/60)*" minutes")

t1=@elapsed dt1=[ensambleMean(meanDist,n1,k,β,-β) for β=0.8:0.1/intervals:0.9]*3/(2*n1)
writedlm( "graphDiams2_50.csv",  dt1, ',')
println("Graph Diameters 50 took "*string(t1)*"s")

t2=@elapsed dt2=[ensambleMean(meanDist,n2,k,β,-β) for β=0.8:0.1/intervals:0.9]*3/(2*n2)
writedlm( "graphDiams2_100.csv",  dt2, ',')
println("Graph Diameters 100 took "*string(t2)*"s")

t3=@elapsed dt3=[ensambleMean(meanDist,n3,k,β,-β) for β=0.8:0.1/intervals:0.9]*3/(2*n3)
writedlm( "graphDiams2_200.csv",  dt3, ',')
println("Graph Diameters 200 took "*string(t3)*"s")

t4=@elapsed dt4=[ensambleMean(meanDist,n4,k,β,-β) for β=0.8:0.1/intervals:0.9]*3/(2*n4)
writedlm( "graphDiams2_400.csv",  dt4, ',')
println("Graph Diameters 400 took "*string(t4)*"s")

t5=@elapsed dt5=[ensambleMean(meanDist,n5,k,β,-β) for β=0.8:0.1/intervals:0.9]*3/(2*n5)
writedlm( "graphDiams2_800.csv",  dt5, ',')
println("Graph Diameters 800 took "*string(t5)*"s")

println("Total computation time: "*string((t1+t2+t3+t4+t5)/60)*" minutes")
