using StaticArrays, Plots, Printf, WriteVTK

using WaterLily, ParametricBodies
using BiotSavartBCs

using CUDA

function main()
    # +++ Simulation parameters +++
    mem = CuArray
    T = Float32
    zeroT = zero(T)
    oneT  = one(T)

    # +++ Controlling flow parameters +++
    N = 96  # number of grid
    g = 1
    Re = 1e5

    # +++ Derived flow parameters +++
    R = N÷4
    U = T(sqrt(g*R))
    ν = U*R/Re
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
            abs(x[3]-center[3]) - R, 
            sqrt((x[1]-center[1])^2 + (x[2]-center[2])^2) - R
        )
    )
    bullet = sphere + cylinder

    body = prolate
    sim = BiotSimulation((N÷2,N÷2,2N), uBC, R; ν, body, mem, T, U)

    vtk_v(a::AbstractSimulation) = a.flow.u/a.U |> Array
    vtk_p(a::AbstractSimulation) = a.flow.p/(0.5a.U^2) |> Array
    vtk_λ₂(a::AbstractSimulation)  = (@inside a.flow.σ[I] = WaterLily.λ₂(I,a.flow.u)*a.L/a.U; a.flow.σ |> Array)
    vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array)
    custom_write_attributes = Dict("u" => vtk_v, "p" => vtk_p, "λ₂" => vtk_λ₂, "d" => vtk_body)

    wr = vtkWriter("FallingBody"; attrib=custom_write_attributes)

    @time for tᵢ in range(0.,20;step=0.1)
        sim_step!(sim,tᵢ;remeasure=true)
        @printf("tU/L= %5.2f\n", sim_time(sim))
        save!(wr, sim)
    end
    close(wr)

end

main()
