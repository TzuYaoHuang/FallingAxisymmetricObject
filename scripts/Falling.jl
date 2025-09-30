using WaterLily, StaticArrays, Plots, Printf
using WriteVTK
using WaterLily: total_force

let
    # +++ Controlling flow parameters +++
    N = 64
    g = 1
    Re = 1e5

    # +++ Derived flow parameters +++
    R = N÷4
    U = sqrt(g*R)
    ν = U*R/Re
    center= SA[0,0,N/3]

    body = AutoBody((x,t) -> √sum(abs2, x.-center)-R)
    sim = Simulation((N÷2,N÷2,2N), (0,0,U), R; ν, body)

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