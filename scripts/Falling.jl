using StaticArrays, Plots, Printf, WriteVTK, CSV, Tables, DataFrames

using WaterLily, ParametricBodies
using BiotSavartBCs

using CUDA

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
    ReList = [9780, 19005, 12474, 101769]
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
    # Read the file
    f = open("../figure/Cd_N$(N).csv")
    bodyName = split(readline(f), ",")[2:end] # read only first line
    seekstart(f) # reset file pointer to beginning
    AllList = Matrix(CSV.read(f, DataFrame, header=false, skipto=2)) # skip header
    simTime = AllList[:,1]
    CdList = AllList[:,2:end]
    close(f)

    # +++ Post-processing
    for I∈CartesianIndices(CdList) CdList[I] = ifelse(abs(CdList[I])>10, NaN, CdList[I]) end
    plot(size=(600,500))
    for (iBody, nameB)∈enumerate(bodyName)
        plot!(simTime, CdList[:,iBody],label=nameB)
    end
    plot!(ylimit=(0,1), xlimit=(simTime[1],simTime[end]))
    plot!(xlabel="rU/R", ylabel="Cd")
    savefig("../figure/CD_N$(N).png")
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
