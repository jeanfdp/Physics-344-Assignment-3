using Combinatorics, Graphs
using StatsBase, Random, Distributions
using DelimitedFiles, GraphIO
using Memoization

@memoize function splitEdges(n)#memoized since can be called multiple times for same size faster
    grd=Graphs.grid([n,n])#generate underlyin grid graph
    gre=unique(map(x->sort(x),collect(Iterators.flatten(map(x->collect(powerset(x, 2, 2)),filter(x->length(x)==4,simplecycles_limited_length(grd,4,4*n^2)))))))#Find all shortcuts by taking all possible subsets with length 2 of length 4 cycles
    return [filter(x->has_edge(grd,x[1],x[2]),gre),filter(x->!has_edge(grd,x[1],x[2]),gre)]#Split up the shortcuts according to their strain, by checking if they are in the underlying graph
end

function ensambleMean(func,n,k,β,μ)
    out=Threads.Atomic{Float64}(0)#Threadsafe total
    edges=splitEdges(n)#Calculate all possible shiortcuts
    (p1,p2)=(1-1/(1+exp(β*(μ))),1-1/(1+exp(β*(μ-1))))#calculate the probability of given edge of strain 0 and 1 existing
    Threads.@threads for i=1:k
        (n1,n2)=(rand(Binomial(length(edges[1]),p1)),rand(Binomial(length(edges[2]),p2)))#Sample binomial distributions to get typical amount of edges for given ensemble
        if n1+n2>0#Then if the graph isnt empty, we apply the supllied function to the graph we can construct from the chosen edges.
            Threads.atomic_add!(out,func(Graphs.SimpleGraphs.SimpleGraphFromIterator(map(x->Graphs.SimpleEdge(x[1],x[2]),vcat(StatsBase.sample(edges[1],n1,replace=false),StatsBase.sample(edges[2],n2,replace=false))))))
        end
    end
    return out[]/(k)#Return mean
end


function binarySearchPercolation(func,n1,n2,k,steps,β1,β2,μ)#Do standard binary search
    for i=1:steps
        Δ=ensambleMean(func,n2,k,(β1+β2)/2,μ)/n2^2-ensambleMean(func,n1,k,(β1+β2)/2,μ)/n1^2
        if Δ>=0
            β2=(β1+β2)/2
        else
            β1=(β1+β2)/2
        end
    end
    return (β1+β2)/2
end

######################
###Bond Percolation###
######################
function maxCompSize(g)
    return float(maximum(map(x->length(x),connected_components(g))))#return amount of vertices in largest component
end


##############################
###Crystal Site Percolation###
##############################
function max4CompSize(g)
    verts=findall(x->x==4,degree(g))#Find all vertices of degree 4
    return maxCompSize(induced_subgraph(g,verts)[1])#Find largest component size in subgraph of vertices with degree 4
end


#################################
###Degenerate site Percolation###
#################################
function max8CompSize(g)
    verts=findall(x->x==8,degree(g))#Find all vertices of degree 8
    return maxCompSize(induced_subgraph(g,verts)[1])#Find largest component size in subgraph of vertices with degree 8
end


#############
###Running###
#############
const n1=50#smaller graph size
const n2=100#larger graph size
const k=1024#amount of graphs used to estimate ensemble mean
const steps=8#Steps in binary search
const interval=0.05#resolution on mu axis



t1=@elapsed dt1=[binarySearchPercolation(maxCompSize,n1,n2,k,steps,20/(1-μ),0,μ) for μ=-3:interval:-interval]
writedlm("negBondPercolation.csv",  dt1, ',')
println("Part 1 of bond percolation took "*string(t1)*"s")

t2=@elapsed dt2=[binarySearchPercolation(maxCompSize,n1,n2,k,steps,-20/μ,0,μ) for μ=1+interval:interval:3]
writedlm("BondPercolation.csv",  dt2, ',')
println("Part 2 of bond percolation took "*string(t2)*"s")


t3=@elapsed dt3=[binarySearchPercolation(max8CompSize,n1,n2,k,steps,0,-20/(1-μ),μ) for μ=-3:interval+interval:-interval]
writedlm("neg8SitePercolation.csv",  dt3, ',')
println("Part 1 of 8-site percolation took "*string(t3)*"s")

t4=@elapsed dt4=[binarySearchPercolation(max8CompSize,n1,n2,k,steps,1,20/μ,μ) for μ=1+interval:interval:3]
writedlm("8SitePercolation.csv",  dt4, ',')
println("Part 2 of 8-site percolation took "*string(t4)*"s")


t5=@elapsed dt5=[binarySearchPercolation(max4CompSize,n1,n2,k,steps,-5,-20,μ) for μ=interval:interval:1-interval]
writedlm("neg4SitePercolation.csv",  dt5, ',')
println("Part 1 of 4-site percolation took "*string(t5)*"s")

t6=@elapsed dt6=[binarySearchPercolation(max4CompSize,n1,n2,k,steps,5,20,μ) for μ=interval:interval:1-interval]
writedlm("4SitePercolation.csv",  dt6, ',')
println("Part 2 of 4-site percolation took "*string(t6)*"s")

t7=@elapsed dt7=[binarySearchPlasma(n1,n2,k,steps,-0.8-1.5/μ,0.8-1.5/μ,μ) for μ=-3:interval:-interval]
writedlm("neg0SitePercolation.csv",  dt7, ',')
println("Part 1 of 0-site percolation took "*string(t7)*"s")

t8=@elapsed dt8=[binarySearchPlasma(n1,n2,k,steps,0.8-1.5/(μ-1),-0.8-1.5/(μ-1),μ) for μ=1+interval:interval:3]
writedlm("0SitePercolation.csv",  dt8, ',')
println("Part 2 of 0-site percolation took "*string(t8)*"s")


println("Total computation time: "*string((t1+t2+t3+t4+t5+t6+t7+t8)/60)*" minutes")
