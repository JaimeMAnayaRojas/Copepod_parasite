using Agents, Random
using InteractiveDynamics
using CairoMakie


mutable struct GrazerCopepod <: AbstractAgent
    id::Int
    pos::Dims{2}
    type::Symbol # :grazer or :copepod
    energy::Float64
    reproduction_prob::Float64
    Δenergy::Float64
end

Grazer(id, pos, energy, repr, Δe) = GrazerCopepod(id, pos, :grazer, energy, repr, Δe)
Copepod(id, pos, energy, repr, Δe) = GrazerCopepod(id, pos, :copepod, energy, repr, Δe)

function initialize_model(;
    n_grazer = 100,
    n_copepods = 50,
    dims = (20, 20), # 2D
    #dims = (100., 100., 50.), # 3D
    regrowth_time = 30,
    Δenergy_grazer = 4,
    Δenergy_copepod = 20,
    grazer_reproduce = 0.04,
    copepod_reproduce = 0.05,
    seed = 23182,
)

    rng = MersenneTwister(seed)
     space = GridSpace(dims, periodic = false) # 
    # Note that the dimensions of the space do not have to correspond to the dimensions
    # of the pathfinder. Discretisation is handled by the pathfinding methods
    #space = ContinuousSpace(dims; periodic = false)
    
    
    
    
    # Model properties contain the phytoplankton as two arrays: whether it is fully grown
    # and the time to regrow. Also have static parameter `regrowth_time`.
    # Notice how the properties are a `NamedTuple` to ensure type stability.
    properties = (
        fully_grown = falses(dims),
        countdown = zeros(Int, dims),
        regrowth_time = regrowth_time,
    )
    model = ABM(GrazerCopepod, space; properties, rng, scheduler = Schedulers.randomly)
    id = 0
    for _ in 1:n_grazer
        id += 1
        energy = rand(1:(Δenergy_grazer*2)) - 1
        grazer = Grazer(id, (rand(0:20), rand(0:20)), energy, grazer_reproduce, Δenergy_grazer)
        add_agent!(grazer, model)
    end
    for _ in 1:n_copepods
        id += 1
        energy = rand(1:(Δenergy_copepod*2)) - 1
        copepod = Copepod(id, (rand(0:20), rand(0:20)), energy, copepod_reproduce, Δenergy_copepod)
        add_agent!(copepod, model)
    end
    for p in positions(model) # random phytoplankton initial growth
        fully_grown = rand(model.rng, Bool)
        countdown = fully_grown ? regrowth_time : rand(model.rng, 1:regrowth_time) - 1
        model.countdown[p...] = countdown
        model.fully_grown[p...] = fully_grown
    end
    return model
end


function grazercopepod_step!(agent::GrazerCopepod, model)
    if agent.type == :grazer
        grazer_step!(agent, model)
    else # then `agent.type == :copepod`
        copepod_step!(agent, model)
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

function copepod_step!(copepod, model)
    walk!(copepod, rand, model)
    copepod.energy -= 1
    agents = collect(agents_in_position(copepod.pos, model))
    dinner = filter!(x -> x.type == :grazer, agents)
    copepod_eat!(copepod, dinner, model)
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

function copepod_eat!(copepod, grazer, model)
    if !isempty(grazer)
        dinner = rand(model.rng, grazer)
        kill_agent!(dinner, model)
        copepod.energy += copepod.Δenergy 
    end
end

function reproduce!(agent, model)
    agent.energy /= 2

    if agent.type == :grazer

        for i in 1:rand(1:10)
            id = nextid(model)
            offspring = GrazerCopepod(
                id,
                agent.pos,
                agent.type,
                agent.energy,
                agent.reproduction_prob,
                agent.Δenergy,
            )
            add_agent_pos!(offspring, model)
        
        end
    else
        for i in 1:rand(1:5)
            id = nextid(model)
            offspring = GrazerCopepod(
                id,
                agent.pos,
                agent.type,
                agent.energy,
                agent.reproduction_prob,
                agent.Δenergy,
            )
            add_agent_pos!(offspring, model)
        
        end
    end

   
    return
end

function phytoplankton_step!(model)
    @inbounds for p in positions(model) # we don't have to enable bound checking
        if !(model.fully_grown[p...])
            if model.countdown[p...] ≤ 0
                model.fully_grown[p...] = true
                model.countdown[p...] = model.regrowth_time
            else
                model.countdown[p...] -= 1
            end
        end
    end
end

model = initialize_model()

offset(a) = a.type == :grazer ? (-0.7, -0.5) : (-0.3, -0.5)
ashape(a) = a.type == :grazer ? :circle : :utriangle
acolor(a) = a.type == :grazer ? RGBAf0(1.0, 1.0, 1.0, 0.8) : RGBAf0(0.2, 0.2, 0.2, 0.8)

phytoplanktoncolor(model) = model.countdown ./ model.regrowth_time

heatkwargs = (colormap = [:brown, :green], colorrange = (0, 1))

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
copepods(a) = a.type == :copepod
count_phytoplankton(model) = count(model.fully_grown)

model = initialize_model()
n = 500
adata = [(grazer, count), (copepods, count)]
mdata = [count_phytoplankton]
adf, mdf = run!(model, grazercopepod_step!, phytoplankton_step!, n; adata, mdata)


function plot_population_timeseries(adf, mdf)
    figure = Figure(resolution = (600, 400))
    ax = figure[1, 1] = Axis(figure; xlabel = "Step", ylabel = "Population")
    grazerl = lines!(ax, adf.step, adf.count_grazer, color = :blue)
    copepodl = lines!(ax, adf.step, adf.count_copepods, color = :orange)
    phytol = lines!(ax, mdf.step, mdf.count_phytoplankton, color = :green)
    figure[1, 2] = Legend(figure, [grazerl, copepodl, phytol], ["Grazers", "Copepods", "Algae"])
    figure
end

plot_population_timeseries(adf, mdf)


n = 2000

model = initialize_model(
    n_grazer = 500,
    n_copepods = 500,
    dims = (500, 500),
    regrowth_time = 30,
    Δenergy_grazer = 4,
    Δenergy_copepod = 20,
    grazer_reproduce = 0.08,
    copepod_reproduce = 0.04,
    seed = 23182,
)


adf, mdf = run!(model, grazercopepod_step!, phytoplankton_step!, n; adata, mdata)

plot_population_timeseries(adf, mdf)


plot(mdf.count_phytoplankton, adf.count_grazer)
plot(adf.count_grazer, adf.count_copepods)
plot(mdf.count_phytoplankton, adf.count_copepods)


figure = Figure(resolution = (600, 400))
ax = figure[1, 1] = Axis(figure; xlabel = "Prey", ylabel = "Predator")
grazerl = lines!(ax, mdf.count_phytoplankton, adf.count_grazer, color = :blue)
copepodl = lines!(ax, adf.count_grazer, adf.count_copepods, color = :orange)
phytol = lines!(ax, mdf.count_phytoplankton, adf.count_copepods,  color = :green)
figure[1, 2] = Legend(figure, [grazerl, copepodl, phytol], ["Grazers vs Phyto", "Copepods vs. Grazer", "Copepod vs. Algae"])
figure


####
