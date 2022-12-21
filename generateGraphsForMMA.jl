using Combinatorics, Graphs
using StatsBase, Random, Distributions
using DelimitedFiles, GraphIO
using Memoization
using MathLink

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

weval(W`MeanGraphDistance2[graph_] := N@Mean[GraphDistance[graph, #1, #2] & @@@RandomChoice[VertexList[graph], {500, 2}]]`)

function biggestDiam(g)
    nm=string(rand(Int))*".el"
    savegraph("./tmp/"*nm,g,EdgeListFormat())
    cmd=W`MeanGraphDistance2[Subgraph[#, First@ConnectedComponents[#]]]&[Graph[Import[FileNameJoin[{Directory[], "tmp", name}]]]]`
    out=weval(cmd,name=nm)
    rm("./tmp/"*nm)
    return float(out)
end

const n1=50
const n2=100
const n3=200
const n4=400
const n5=800
const n6=1200
const k=128
const intervals=10

println("Starting!")

# t1=@elapsed dt1=[ensambleMean(biggestDiam,n1,k,β,-β) for β=0.5:0.5/intervals:1]*3/(2*n1)
# writedlm( "graphDiams50.csv",  dt1, ',')
# println("Graph Diameters 50 took "*string(t1)*"s")

# t2=@elapsed dt2=[ensambleMean(biggestDiam,n2,k,β,-β) for β=0.5:0.5/intervals:1]*3/(2*n2)
# writedlm( "graphDiams100.csv",  dt2, ',')
# println("Graph Diameters 100 took "*string(t2)*"s")

# t3=@elapsed dt3=[ensambleMean(biggestDiam,n3,k,β,-β) for β=0.5:0.5/intervals:1]*3/(2*n3)
# writedlm( "graphDiams200.csv",  dt3, ',')
# println("Graph Diameters 200 took "*string(t3)*"s")

# t4=@elapsed dt4=[ensambleMean(biggestDiam,n4,k,β,-β) for β=0.5:0.5/intervals:1]*3/(2*n4)
# writedlm( "graphDiams400.csv",  dt4, ',')
# println("Graph Diameters 400 took "*string(t4)*"s")

t5=@elapsed dt5=[ensambleMean(biggestDiam,n5,k,β,-β) for β=0.5:0.5/intervals:1]*3/(2*n5)
writedlm( "graphDiams800.csv",  dt5, ',')
println("Graph Diameters 800 took "*string(t5)*"s")

t6=@elapsed dt6=[ensambleMean(biggestDiam,n6,k,β,-β) for β=0.5:0.5/intervals:1]*3/(2*n6)
writedlm( "graphDiams1200.csv",  dt6, ',')
println("Graph Diameters 1200 took "*string(t6)*"s")

println("Total computation time: "*string((t5+t6)/60)*" minutes")
# println("Total computation time: "*string((t1+t2+t3+t4+t5)/60)*" minutes")