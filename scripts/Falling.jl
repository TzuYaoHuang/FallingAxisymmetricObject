using StaticArrays, Plots, Printf, WriteVTK, CSV, Tables, DataFrames
using Statistics

using WaterLily, ParametricBodies
using BiotSavartBCs

# using CUDA

# Include the mirror in x and y plane
import BiotSavartBCs: interaction,image,symmetry
@inline function symmetry(ω,T,args...) # overwrite to add image influences
    T₁,sgn₁ = image(T,size(ω),-1)
    T₂,sgn₂ = image(T,size(ω),-2)
    T₁₂,_   = image(T₁,size(ω),-2)
    return interaction(ω,T,args...)+sgn₁*interaction(ω,T₁,args...)+
        sgn₂*(interaction(ω,T₂,args...)+sgn₁*interaction(ω,T₁₂,args...))
end

function determineMode(args)
    if length(args) ==0 
        error("Please specify the mode!")
    end
    modeFlag = args[1]
    simulation = modeFlag=="run" || modeFlag=="all"
    postProcess= modeFlag=="pp" || modeFlag=="all"
    return simulation,postProcess
end

function run(N; mem=CuArray, T=Float32, tEnd=30)
    # +++ Simulation parameters +++
    zeroT = zero(T)
    oneT  = one(T)

    # +++ Controlling flow parameters +++
    U = T(1)

    # +++ Derived flow parameters +++
    R = N÷4
    D = 2R
    A = π*R^2/4
    center= SA{T}[0,0,2N/3]
    uBC(i,x,t) = ifelse(i==3, U, zeroT)

    # +++ Body Definition -- Sphere
    sphere = AutoBody((x,t) -> √sum(abs2, x.-center)-R)
    # +++ Body Definition -- Prolate
    γ = T(1.5)
    function prolateSDF(x,t)
        X,Y,Z = @. x-center
        F = (X/R)^2 + (Y/R)^2 + (Z/(γ*R))^2 - oneT
        gradnorm = sqrt((2X/R^2)^2 + (2Y/R^2)^2 + (2Z/(γ*R)^2)^2)
        return F/gradnorm
    end
    prolate  = AutoBody(prolateSDF)
    # +++ Body Definition -- Bullet
    cylinder = AutoBody(
        (x,t) -> max(
            abs(x[3]-(center[3]+R)) - R, 
            sqrt((x[1]-center[1])^2 + (x[2]-center[2])^2) - R
        )
    )
    bullet = sphere + cylinder
    # +++ Body Definition -- NACA teardrop
    NACA_twodigits = 31.0559
    thickness = T(NACA_twodigits/100)
    maxThick_x = T(0.3)
    NACA(s) = thickness*5*(0.2969f0s-0.126f0s^2-0.3516f0s^4+0.2843f0s^6-0.1036f0s^8)
    curve(s,t) = 2R/thickness*SA[(1-s)^2-maxThick_x,NACA(1-s)] .+ SA[center[3],0]
    revolve(x::SVector{3},t) = SA[x[3],hypot(x[1],x[2])] # revolve around x[3]-axis
    teardrop = ParametricBody(curve,HashedLocator(curve,(0,1);T,mem);map=revolve,ndims=3)

    # +++ List of all body +++
    bodies = [sphere, prolate, bullet, teardrop]
    ReList = [9994, 12241, 8316, 17311]
    νList = U*D./ReList
    bodyName=["Sphere", "Ellipsoid", "Bullet", "Teardrop"]
    NBody = length(bodies)

    # ++++ VTK Writer function
    vtk_v(a::AbstractSimulation) = a.flow.u/a.U |> Array
    vtk_p(a::AbstractSimulation) = a.flow.p/(0.5a.U^2) |> Array
    vtk_λ₂(a::AbstractSimulation)  = (@inside a.flow.σ[I] = WaterLily.λ₂(I,a.flow.u)*a.L/a.U; a.flow.σ |> Array)
    vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array)
    custom_write_attributes = Dict("u" => vtk_v, "p" => vtk_p, "λ₂" => vtk_λ₂, "d" => vtk_body)

    # ++++ Generate Simulation
    simTime = 0:0.1:tEnd; NTime = length(simTime)
    CdList = zeros(NTime,NBody) 
    for (iBody, body)∈enumerate(bodies)
        println("$(bodyName[iBody]) is falling now.")

        # disable BiotSavartBCs in -x, and -y faces
        sim = BiotSimulation((N÷2,N÷2,3N), uBC, R; ν=νList[iBody], body, mem, T, U, nonbiotfaces=(-1,-2))
        wr = vtkWriter("Falling_$(bodyName[iBody])_N$(N)"; attrib=custom_write_attributes)

        # Running Simulation!
        for (iTime,tᵢ) in enumerate(simTime)
            sim_step!(sim,tᵢ;remeasure=true)
            @printf("tU/L= %5.2f\n", sim_time(sim)); flush(stdout)
            save!(wr, sim)
            CdList[iTime,iBody] = -WaterLily.total_force(sim)[3]/(0.5*sim.U^2*A)
        end
        close(wr)
    end
    CSV.write("../figure/Cd_N$(N).csv",  Tables.table(hcat(collect(simTime), CdList)), header=["t", bodyName...])
end

function postProcess(N)
    # Read and parse data
    df = CSV.read("../figure/Cd_N$(N).csv", DataFrame)
    bodyNames = names(df)[2:end]  # All columns except time
    simTime = df.t
    CdData = Matrix(df[:, 2:end])
    
    # Data cleaning and validation
    for I ∈ CartesianIndices(CdData) 
        CdData[I] = ifelse(abs(CdData[I]) > 10 || isnan(CdData[I]), NaN, CdData[I]) 
    end
    
    # Scientific analysis
    println("\n" * "="^80)
    println("DRAG COEFFICIENT ANALYSIS - Resolution N=$N")
    println("="^80)
    
    # Terminal values analysis (last 20% of simulation)
    terminal_idx = max(1, round(Int, 0.8 * length(simTime))):length(simTime)
    terminal_stats = DataFrame(
        Body = String[],
        Mean_Cd = Float64[],
        Std_Cd = Float64[],
        Min_Cd = Float64[],
        Max_Cd = Float64[],
        Terminal_Cd = Float64[]
    )
    
    for (i, bodyName) in enumerate(bodyNames)
        terminal_data = filter(!isnan, CdData[terminal_idx, i])
        if !isempty(terminal_data)
            push!(terminal_stats, (
                bodyName,
                round(mean(terminal_data), digits=4),
                round(std(terminal_data), digits=4),
                round(minimum(terminal_data), digits=4),
                round(maximum(terminal_data), digits=4),
                round(CdData[end, i], digits=4)
            ))
        end
    end
    
    println(terminal_stats)
    
    # Create professional visualization
    # Color palette - scientific and colorblind-friendly
    colors = [:steelblue, :darkorange, :forestgreen, :crimson, :purple, :brown]
    line_styles = [:solid, :dash, :dot, :dashdot]
    
    # Main time series plot
    p1 = plot(
        size = (900, 600),
        dpi = 300,
        background_color = :white,
        foreground_color = :black,
        grid = true,
        gridwidth = 1,
        gridcolor = :lightgray,
        gridstyle = :dot,
        framestyle = :box,
        thickness_scaling = 1.2
    )
    
    for (i, bodyName) in enumerate(bodyNames)
        valid_idx = .!isnan.(CdData[:, i])
        if any(valid_idx)
            plot!(p1, simTime[valid_idx], CdData[valid_idx, i],
                label = bodyName,
                color = colors[mod1(i, length(colors))],
                linestyle = line_styles[mod1(i, length(line_styles))],
                linewidth = 2.5,
                markershape = :none
            )
        end
    end
    
    plot!(p1,
        xlabel = "Dimensionless Time (tU/R)",
        ylabel = "Drag Coefficient (C_d)",
        title = "Drag Coefficient Evolution - N=$N Grid Resolution",
        titlefontsize = 14,
        guidefontsize = 12,
        tickfontsize = 10,
        legendfontsize = 10,
        legend = :topright,
        xlims = (simTime[1], simTime[end]),
        ylims = (0, maximum(filter(!isnan, CdData)) * 1.1),
        left_margin = 5Plots.mm,
        bottom_margin = 5Plots.mm,
        top_margin = 3Plots.mm,
        right_margin = 3Plots.mm
    )
    
    # Add horizontal lines for terminal values
    for (i, row) in enumerate(eachrow(terminal_stats))
        hline!(p1, [row.Mean_Cd], 
            color = colors[mod1(i, length(colors))],
            linestyle = :dashdotdot,
            alpha = 0.6,
            linewidth = 1,
            label = ""
        )
    end
    
    # Convergence analysis plot
    p2 = plot(
        size = (900, 400),
        dpi = 300,
        background_color = :white,
        grid = true,
        gridwidth = 1,
        gridcolor = :lightgray,
        gridstyle = :dot,
        framestyle = :box
    )
    
    # Moving average for convergence
    window_size = max(5, length(simTime) ÷ 20)
    for (i, bodyName) in enumerate(bodyNames)
        valid_data = CdData[:, i]
        valid_idx = .!isnan.(valid_data)
        if sum(valid_idx) > window_size
            moving_avg = [mean(valid_data[max(1,j-window_size):j]) 
                         for j in window_size:length(valid_data) if valid_idx[j]]
            time_avg = simTime[window_size:end][valid_idx[window_size:end]]
            
            plot!(p2, time_avg, moving_avg,
                label = "$bodyName (Moving Avg)",
                color = colors[mod1(i, length(colors))],
                linewidth = 2,
                alpha = 0.8
            )
        end
    end
    
    plot!(p2,
        xlabel = "Dimensionless Time (tU/R)",
        ylabel = "Cd (Moving Average)",
        title = "Convergence Analysis - Moving Average",
        titlefontsize = 12,
        guidefontsize = 11,
        tickfontsize = 9,
        legendfontsize = 9,
        legend = :right,
        xlims = (simTime[1], simTime[end])
    )
    
    # Combined layout
    final_plot = plot(p1, p2, 
        layout = (2, 1), 
        size = (900, 1000),
        plot_title = "Computational Fluid Dynamics Analysis",
        plot_titlefontsize = 16
    )
    
    # Save in multiple formats
    savefig(final_plot, "../figure/CD_N$(N)_analysis.png")
    
    # Export terminal statistics
    CSV.write("../figure/terminal_stats_N$(N).csv", terminal_stats)
    
    println("\nAnalysis complete! Files saved:")
    println("  - Figures: CD_N$(N)_analysis.png")
    println("  - Statistics: terminal_stats_N$(N).csv")
    println("="^80)
    
    return final_plot, terminal_stats
end

function main()
    simu,pp = determineMode(ARGS)
    N = parse(Int, ARGS[2])

    simu && run(N)
    pp && postProcess(N)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
