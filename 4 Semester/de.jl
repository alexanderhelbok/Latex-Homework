using DataFrames, Statistics
using StatsBase
using GLMakie

function chisq(obs, exp; sigma=nothing, dof=nothing, pcount=nothing)
    """
    Calculate chi-squared value for a given set of observed and expected values.
    Parameters:
    obs: observed values
    exp: expected values
    sigma: optional error on observed values (default nothing)
    dof: degrees of freedom (default 0)
    pcount: number of parameters (default 0)
    Returns:
    chi-squared value
    """
    if dof === nothing && pcount === nothing
        dof = 1
    elseif dof === nothing
        dof = length(obs) - pcount
    end

    if sigma === nothing
        return sum((obs - exp).^2)/dof
    else
        # replace 0 values in sigma with 1
        sigma[sigma .== 0] .= 1
        return sum(((obs - exp)./sigma).^2)/dof
    end
end

# PythonPlot.matplotlib.style.use("test.mplstyle")

# include("Source.jl")

@. model(x, p) = p[1] * sin(p[2] * x + p[3]) + p[4]

xdata = range(0, 10, 20)
ydata = model(xdata, [1, 1, -0.3, 0.4]) + 0.05 * randn(length(xdata))


function de(fit, xdata, ydata, bounds; mut=0.8, crossp=0.7, popsize=20, its=1000, fobj=chisq, sigma = nothing, seed=nothing)
    # initial checks
    if length(xdata) != length(ydata)
        error("xdata and ydata must have the same length")
    end
    # set seed for reproducibility
    if seed !== nothing
        Random.seed!(seed)
    end
    

    dimensions = length(bounds)
    # create population with random parameters (between 0 and 1)
    pop = rand(popsize, dimensions)
    # scale parameters to the given bounds
    bounds = permutedims(hcat(bounds...))
    min_b, max_b = bounds[:, 1], bounds[:, 2]
    diff = abs.(max_b - min_b)
    pop_denorm = min_b' .+ pop .* diff'
    # calculate fitness (higher is worse)
    fitness = [fobj(ydata, fit(xdata, p), sigma=sigma) for p in eachrow(pop_denorm)]
    # sort by fitness and get best (lowest) one
    best_idx = argmin(fitness) 
    best = pop_denorm[best_idx, :]
    # start evolution

    pop_t = reshape(zeros(popsize * dimensions * its), popsize, dimensions, its)
    pop_t[:, :, 1] = pop_denorm[:, :]
    for i in 1:its
        if i != 1
            pop_t[:, :, i] = pop_t[:, :, i-1]
        end
        for j in 1:popsize
        # select three random vector index positions (not equal to j) without replacement
            idxs = filter(x -> x != j, 1:popsize)
            mu = @view pop[idxs[sample(1:popsize-1, 3, replace=false)], :]
            # create a mutant by adding random scaled difference vectors
            mutant = @. mu[1, :] + mut * (mu[2, :] - mu[3, :])
            mutant = clamp.(mutant, 0, 1)
            # randomly create a crossover mask
            cross_points = rand(dimensions) .< crossp
            # construct trial vector by mixing the mutant and the current vector
            trial = [if cross_points[k] 
                        mutant[k] 
                    else 
                        pop[j, k] 
                    end for k in 1:dimensions]

            trial_denorm = min_b .+ trial .* diff

            # calculate fitness
            f = fobj(ydata, fit(xdata, trial_denorm), sigma=sigma)
            # replace the current vector if the trial vector is better
            if f < fitness[j]
                fitness[j] = f
                pop[j, :] = trial
                pop_t[j, :, i] = trial_denorm
                if f < fitness[best_idx]
                    best_idx = j
                    best = trial_denorm
                end
            end
        end
    end
    return pop_t, fitness
    # return best, fitness[best_idx]
    # return min_b .+ pop .* diff, fitness
end

best, fitness = de(model, xdata, ydata, [[-2, 2], [-2, 2], [-2, 2], [-2, 2]], mut=0.8, its=1000, popsize=10)

iterator = Observable(1)

nframes = 1000


x = [0:0.1:10;]
y = @lift model(x, best[1, :, $iterator])

begin
fig, ax, s1 = scatter(xdata, ydata, color = :red)
lines!(x, model(x, [1, 1, -0.3, 0.4]), linewidth = 10, color = :magenta, alpha = 0.5, linestyle = :dash)
for i in 1:10
    lines!(x, @lift(model(x, best[i, :, $iterator])), linewidth = 1)
end

record(fig, "time_animation.mp4", 1:nframes;
        framerate = 100) do t
    iterator[] = t
end 
end

