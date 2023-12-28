using Measurements, CSV, DataFrames

function nom(x)
    return Measurements.value(x)
end

function str2meas(s)
    s = replace(s, "+/-" => "±")
    # check for nan
    if occursin("nan", s)
        nom = parse(Float64, split(s, "±")[1])
        return measurement(nom, 0)
    end
    result = parse(Measurement{Float64}, s)
    return result
end

function load_data(run)
    # load data
    data = DataFrame()
    for file in readdir("data/RUN$run")
        tempdf = CSV.read("data/RUN$run/$file", DataFrame)
    
        tempdf.Sample = [file[1:end-4] for i in 1:size(tempdf, 1)]
        # convert columns to Measurement
        tempdf.power = nom.(str2meas.(tempdf.power))
        tempdf.Qc = nom.(str2meas.(tempdf.Qc))
        tempdf.Qint = nom.(str2meas.(tempdf.Qint))
        tempdf.num_phot = nom.(str2meas.(tempdf.num_phot))
        # duplicate Res
        tempdf.Resint = tempdf.Res
        # convert Res to string
        tempdf.Res = string.(tempdf.Res)
        # replace Res 
        tempdf.Res = replace!(tempdf.Res, "1" => "C", "2" => "D2", "3" => "D1")
        sort!(tempdf)
    
        data = vcat(data, tempdf)
    end
    return data
end

function newleg(elems::Function, labels, x, y; fig = figure[1, 1], xgap = 0.05, xmarkergap = 0.04, ygap = 0.05, labelsize = 22)
	# create a temporary axis (necessary for nonlinear axes) to get the relative projection
	tempax = Axis(fig, backgroundcolor = :transparent)
	hidedecorations!(tempax)
	hidespines!(tempax)
	relative_projection = Makie.camrelative(tempax.scene)
	# draw legend markers
	elems(relative_projection, x, y, ygap, xmarkergap)
	# draw legend labels
	for (i, label) in enumerate(labels)
		text!(relative_projection, 
			"$label", 
			position=Point2f(x + xgap, y - (i-1)*ygap), 
			align = (:left, 0.4),
			fontsize = labelsize,
			font = "Times New Roman")
	end
end

function mylegend(figure, elems, labels, x, y; fig = figure[1, 1], rgs...)
	tempax = Axis(fig)
	hidedecorations!(tempax)
	hidespines!(tempax)
	leg_projetion = campixel(tempax.scene)
	@lift translate!(leg_projetion, Vec2f($(figure.scene.camera.resolution)[1]*x, $(figure.scene.camera.resolution)[2]*y))
	Legend(leg_projetion, elems, labels; rgs...)
end

tnr = "Times New Roman"

modes = ["D1", "D2", "C"]
Markers = [:+, :diamond, :circle]
c = [:Blue, :Green, :Red]
# dict for modes with markers and Colors
modedf = DataFrame(mode = modes, marker = Markers, color = c)

modedf

Samples = ["Nb1", "Nb2", "NbTa1", "NbTa2", "Al1", "Al2"]
c = [:darkgreen, :limegreen, :dodgerblue3, :navyblue, :crimson, :brown1]
Sampledf = DataFrame(Sample = Samples, color = c)

Upper(x) = @. 1.0039 - 0.1639*x + 1.8336*x^2 - 2.0191*x^3 - 0.5006*x^4 + 2.5454*x^5 - 2.0145*x^6 + 0.7377*x^7 - 0.1319*x^8 + 0.0093*x^9
Lower(x) = @. 1.0062349088 - 0.2174154056*x + 0.0374311587*x^2 - 0.0039328246*x^3 + 0.000209708*x^4 - 4.6263e-6*x^5 + 3.33e-8*x^6

kwargs = (xminorticksvisible = true,
    yminorticksvisible = true,
    spinewidth = 2,
    xminorticks = IntervalsBetween(4),
    yminorticks = IntervalsBetween(4),
    xtickwidth = 2,
    ytickwidth = 2,
    xticksize = -14,
    yticksize = -14,
    xminorticksize = -7,
    yminorticksize = -7,
    xticksmirrored = true,
    yticksmirrored = true,
    xgridvisible = false,
    ygridvisible = false,
    xgridwidth = 2,
    ygridwidth = 2,
    xticklabelsize = 20,
    yticklabelsize = 20,
    xticklabelfont = "Times New Roman",
    yticklabelfont = "Times New Roman",
    xlabelfont = "Times New Roman",
    xlabelsize = 24,
    ylabelfont = "Times New Roman",
    ylabelsize = 24,
    xlabelpadding = 10,
    ylabelpadding = 10)

legargs = (labelfont = "Times New Roman", 
    labelsize = 20, 
    margin = ones(4).*18,
    bgcolor = :transparent)