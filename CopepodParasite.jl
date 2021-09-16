#TO DO:
# Infected copepods?
# Incoming Parasites? New function? maybe hatch function or another reproduce 
# Influence on infected Copepod behaviour apart from "moving more"
#Line 127, || connection?
# Statistics: add infected copepods somehow 
# only differentiate between phytoplankton and no phytoplankton, count biomass as patcehes where phytocplancton is found
#

using Agents
using Random
pwd()

cd("/home/jaime/Dropbox/Projects_JM/Muenster/Marvin/") # set the working directory
#cd("/home/mao/Dropbox/Projects_JM/Muenster/Marvin/") # set the working directory
readdir() # check the files in the data/ folder


mutable struct CopepodGrazerParasite <: AbstractAgent
    id::Int #Id of the Agent
    pos::Dims{2} #position in the Space
    type::Symbol # :Copepod or :Parasite or :Grazer 
    energy::Float64 
    reproduction_prob::Float64 
    Δenergy::Float64
    infected::Bool
end

function Copepod(id, pos, energy, repr, Δe)
    CopepodGrazerParasite(id, pos, :copepod, energy, repr, Δe, :false)
end

function Parasite(id, pos, energy, repr, Δe)
    CopepodGrazerParasite(id, pos, :parasite, energy, repr, Δe,:false)
end

function Grazer(id, pos, energy, repr, Δe)
    CopepodGrazerParasite(id, pos, :grazer, energy, repr, Δe, :false)
end

function initialize_model(;
    n_copepod = 10,
    n_grazer = 100, # Grazer being Chydoidae, Daphniidae and Sididae 
    n_parasite = 100, #continuous stream of "newly introduced parasites": x amount of bird introduce each: 8000 eggs, only 20% hatch (Merle), check literature 
    dims = (20, 20), #2 Dimensional space
    regrowth_time = 30, #regrowth time of Phytoplankton in "steps"
    Δenergy_copepod = 20, #??? 
    Δenergy_grazer = 20, #??? 
    Δenergy_parasite = 10,#???   
    copepod_reproduce = 0.05, #changes if infected, see copepod_reproduce function
    grazer_reproduce = 0.05, #are not infected -> steady reproduction rate
    parasite_reproduce = 0, 
    seed = 23182,    
)

    rng = MersenneTwister(seed) #MersenneTwister: pseudo random number generator
    space = GridSpace(dims, periodic = false)
    # Model properties contain the phytoplankton as two arrays: whether it is fully grown
    # and the time to regrow. Also have static parameter `regrowth_time`. 
    properties = (
        fully_grown = falses(dims),
        countdown = zeros(Int, dims),
        regrowth_time = regrowth_time,
    )
    model = ABM(CopepodGrazerParasite, space; properties, rng, scheduler = Schedulers.randomly)
    id = 0
    for _ in 1:n_copepod
        id += 1
        energy = rand(1:(Δenergy_copepod*2)) - 1  
        copepod = Copepod(id, (0, 0), energy, copepod_reproduce,Δenergy_copepod)   #all agents start at (0,0)??
        add_agent!(copepod, model)
    end

    for _ in 1:n_grazer
        id += 1
        energy = rand(1:(Δenergy_grazer*2)) - 1  #whats the maximum energy? 
        grazer = Grazer(id, (0,0), energy, grazer_reproduce,Δenergy_grazer)
        add_agent!(grazer, model)
    end
    
    for _ in 1:n_parasite
        id += 1
        energy = rand(1:(Δenergy_parasite*2)) - 1
        parasite = Parasite(id, (0, 0), energy, parasite_reproduce,Δenergy_parasite)
        add_agent!(parasite, model)
    end

    for p in positions(model) #random Phytoplankton initial growth
        fully_grown = rand(model.rng, Bool) #boolean!
        countdown = fully_grown ? regrowth_time : rand(model.rng, 1:regrowth_time) - 1 #if fully grown, use full regrowth timer, else use random number
        model.countdown[p...] = countdown 
        model.fully_grown[p...] = fully_grown #set countdown or fully grown to their respective positions on the map
    end
    return model
end

function copepod_eat!(copepod, food, model)  #copepod eat around their general vicinity 
    if !isempty(food)
        if food == rand(model.rng, food)
            kill_agent!(food, model)
            if food.type == :grazer
                copepod.energy += copepod.Δenergy
            else
                copepod.infected == true
            end
        end
    end
end

function model_step!(agent::CopepodGrazerParasite, model)
    if agent.type == :copepod
        copepod_step!(agent, model)
    elseif agent.type == :grazer 
        grazer_step!(agent, model)
        else 
        parasite_step!(agent, model)
    end
end

function parasite_step!(parasite, model) #in lab: 2 days max (Parasites move really quickly, maybe even follow copepods), copepod 4 days max without food 
    walk!(parasite, rand, model)
    parasite.energy -= 1       
    if parasite.energy < 0
        kill_agent!(parasite, model)
        return
    end
end

function grazer_step!(grazer, model) 
    walk!(grazer, rand, model)
    grazer.energy -= 1 
    grazer_eat!(grazer, model)
    if grazer.energy < 0
        kill_agent!(grazer, model)
        return
    end
    if rand(model.rng) <= grazer.reproduction_prob
        reproduce!(grazer, model)
    end
end

function copepod_step!(copepod, model) #Copepod is able to detect pray at 1mm (parasties want to stay in that vicinity)
    walk!(copepod, rand, model)    #if infected move more around the map 
    if copepod.infected == true
        walk!(copepod,rand, model)
    end 
    copepod.energy -= 1
    agents = collect(agents_in_position(copepod.pos, model)) #does this collect all agents? 
    food = filter!(x -> x.type == :parasite || x.type == :grazer, agents) #add parasites?, maybe && verknüpfung
    #type[findall(x.type .== :parasite, :]
    copepod_eat!(copepod, food, model)  
    if copepod.energy < 0
        kill_agent!(copepod, model)
        return
    end
    if rand(model.rng) <= copepod.reproduction_prob 
        reproduce!(copepod, model)
    end
end

function grazer_eat!(grazer, model)        
    if model.fully_grown[grazer.pos...]
        grazer.energy += grazer.Δenergy
        model.fully_grown[grazer.pos...] = false
    end
end


function reproduce!(agent, model)   #time to grow up 
    if agent.type == :copepod && agent.infected == true 
    else 
    agent.energy /= 2
    id = nextid(model)
    offspring = CopepodGrazerParasite(
        id,
        agent.pos,
        agent.type,
        agent.energy,
        agent.reproduction_prob,
        agent.Δenergy,
        agent.infected
    )
    add_agent_pos!(offspring, model)
    return
    end
end


#checks every position on the model if the phytoplankton is grown fully
function phytoplankton_step!(model)  
    @inbounds for p in positions(model) #@inbounds: loop macro, skip bound checks!
        if !(model.fully_grown[p...])
            if model.countdown[p...] ≤ 0  #time to regrow -> if <0 its fully grown
                model.fully_grown[p...] = true
                model.countdown[p...] = model.regrowth_time #reset timer
            else
                model.countdown[p...] -= 1 #takes one step off time to regrow
            end
        end
    end
end

model = initialize_model() 


using InteractiveDynamics
using CairoMakie


function offset(a)
    a.type == :copepod ? (-0.7, -0.5) : (-0.3, -0.5)
end

function ashape(a)
    if a.type == :copepod 
        :circle 
    elseif a.type == :grazer
        :utriangle
    else
        :hline
    end
end

function acolor(a)
    if a.type == :copepod
        :black 
    elseif (a.type == :copepod) && (a.infected == true)
        :red
    elseif a.type == :grazer 
        :yellow
    else
        :magenta
    end
end

phytoplanktoncolor(model) = model.countdown ./ model.regrowth_time

heatkwargs = (colormap = [:darkseagreen1, :darkgreen], colorrange = (0, 1))

plotkwargs = (
    ac = acolor,
    as = 15,
    am = ashape,
    offset = offset,
    heatarray = phytoplanktoncolor,
    heatkwargs = heatkwargs,
)

fig, _ = abm_plot(model; plotkwargs...)
fig

grazer(a) = a.type == :grazer
copepod(a) = a.type == :copepod
parasite(a) = a.type == :parasite
count_phytoplankton(model) = count(model.fully_grown)


model = initialize_model()

model = initialize_model(;
n_copepod = 100,
n_grazer = 500, # Grazer being Chydoidae, Daphniidae and Sididae 
n_parasite = 100, #continuous stream of "newly introduced parasites": x amount of bird introduce each: 8000 eggs, only 20% hatch (Merle), check literature 
dims = (20, 20), #2 Dimensional space
regrowth_time = 1, #regrowth time of Phytoplankton in "steps"
Δenergy_copepod = 20, #??? 
Δenergy_grazer = 5, #??? 
Δenergy_parasite = 10,#???   
copepod_reproduce = 0.05, #changes if infected, see copepod_reproduce function
grazer_reproduce = 0.005, #are not infected -> steady reproduction rate
parasite_reproduce = 0, 
seed = 23182 
)
n = 500
adata = [(grazer, count), (copepod, count), (parasite, count)]
mdata = [count_phytoplankton]
adf, mdf = run!(model, model_step!, phytoplankton_step!, n; adata, mdata)

adf

mdf

function plot_population_timeseries(adf, mdf)

    figure = Figure(resolution = (600, 400))
    ax = figure[1, 1] = Axis(figure; xlabel = "Step", ylabel = "Population")
    grazer = lines!(ax, adf.step, adf.count_grazer, color = :yellow)
    copepod = lines!(ax, adf.step, adf.count_copepod, color = :black)
    parasite = lines!(ax, adf.step, adf.count_parasite, color = :red)
    phytoplankton = lines!(ax, mdf.step, mdf.count_phytoplankton, color = :green)
    figure[1, 2] = Legend(figure, [grazer, copepod, parasite, phytoplankton], ["Grazers", "Copepods", "Parasites", "Phytoplankton"])
    figure
end


model = initialize_model(;
n_copepod = 500,
n_grazer = 1000, # Grazer being Chydoidae, Daphniidae and Sididae 
n_parasite = 100, #continuous stream of "newly introduced parasites": x amount of bird introduce each: 8000 eggs, only 20% hatch (Merle), check literature 
dims = (20, 20), #2 Dimensional space
regrowth_time = 3, #regrowth time of Phytoplankton in "steps"
Δenergy_copepod = 100, #??? 
Δenergy_grazer = 2, #??? 
Δenergy_parasite = 1,#???   
copepod_reproduce = 0.05, #changes if infected, see copepod_reproduce function
grazer_reproduce = 0.005, #are not infected -> steady reproduction rate
parasite_reproduce = 0, 
seed = 23182 
)
n = 500
adata = [(grazer, count), (copepod, count), (parasite, count)]
mdata = [count_phytoplankton]
adf, mdf = run!(model, model_step!, phytoplankton_step!, n; adata, mdata)

plot_population_timeseries(adf, mdf)


abm_video(
    "copepod.mp4",
    model,
    model_step!, 
    phytoplankton_step!;
    frames = 150,
    framerate = 8,
    plotkwargs...,
)
