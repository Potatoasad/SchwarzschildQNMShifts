module SchwarzschildQNMShifts

using ContourIntegrals
using LeaverSchwQNM

## Define the wavefunction and integrand for Schwarzschild shifts
ψ = RadialMode(2,2,0)
ψ2 = RadialMode(2,2,1)
δV(r,k) = cos(2*π*k/r)

include("CosineShift.jl")

EquationIntegrand = ZeroOrderOperator(ψ,ψ2)
EquationIntegrandSwitched = ZeroOrderOperator(ψ2,ψ)

function norm(integrand)
    δr = 0.4
    bottomright = 1.0+δr-δr*im
    bottomleft = 1.0-δr-δr*im
    c1 = SemiInfiniteLine(bottomright,bottomright+im*1.0,false)
    c2 = LineSegment(bottomright,bottomleft)
    c3 = SemiInfiniteLine(bottomleft,bottomleft+im*1.0,true)

    a1,e1 = Integrate(integrand,c1)
    a2,e2 = Integrate(integrand,c2)
    a3,e3 = Integrate(integrand,c3)
    (a1+a2+a3), sqrt(e1^2 + e2^2 + e3^2)
end

## Do the Integral to check if the norm is self adjoint
ψ(1+0.1+200*im)|> display
ψ2(1+0.1+200*im)|> display
ψ(1-0.1+200*im)|> display
ψ2(1-0.1+200*im)|> display


norm(EquationIntegrand)|> display
norm(EquationIntegrandSwitched)|> display




## Define Domain
δr = 0.4
bottomright = 1.0+δr-δr*im
bottomleft = 1.0-δr-δr*im
c1 = SemiInfiniteLine(bottomright,bottomright+im*1.0,false)
c2 = LineSegment(bottomright,bottomleft)
c1 = SemiInfiniteLine(bottomleft,bottomleft+im*1.0,true)


end # module
