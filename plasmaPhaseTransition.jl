using Combinatorics, Graphs
using StatsBase, Random, Distributions
using DelimitedFiles, GraphIO

function maxCompSize(g)
    return float(maximum(map(x->length(x),connected_components(g))))#return amount of vertices in largest component
end

function emptyComponentSize(n,k,β,μ)#Estimate the ensemble mean of the induced 0-site component size
    out=Threads.Atomic{Float64}(0)
    grd=Graphs.grid([n,n])
    p=(1/(1+exp(β*μ)))^4 *(1/(1+exp(β*(μ-1))))^4#Calculate probability of a vertex having no edges
    Threads.@threads for i=1:k
        m=rand(Binomial(n^2,p))#Sample binomial distribution to get typical amount of edges for given ensemble
        Threads.atomic_add!(out,maxCompSize(induced_subgraph(grd,StatsBase.sample(vertices(grd),m,replace=false))[1]))
    end
    return out[]/(k)
end

function binarySearchPlasma(n1,n2,k,steps,β1,β2,μ)
    for i=1:steps
        Δ=emptyComponentSize(n2,k,(β1+β2)/2,μ)/n2^2-emptyComponentSize(n1,k,(β1+β2)/2,μ)/n1^2
        if Δ>=0
            β2=(β1+β2)/2
        else
            β1=(β1+β2)/2
        end
    end
    return (β1+β2)/2
end



n1=50
n2=100
k=256
steps=8
interval=0.1


t1=@elapsed dt1=[binarySearchPlasma(n1,n2,k,steps,-0.8-1.5/μ,0.8-1.5/μ,μ) for μ=-3:interval:-interval]
writedlm("neg0SitePercolation.csv",  dt1, ',')

t2=@elapsed dt2=[binarySearchPlasma(n1,n2,k,steps,0.8-1.5/(μ-1),-0.8-1.5/(μ-1),μ) for μ=1+interval:interval:3]
writedlm("0SitePercolation.csv",  dt2, ',')

