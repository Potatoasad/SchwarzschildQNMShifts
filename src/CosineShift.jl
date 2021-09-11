function ShiftedOperator(ψs,δV,k)
    dψ = ∂r(ψ);
    d2ψ = ∂r(dψ);
    f = let ψ=ψs, s=ψ.s, l=ψ.l, ρ = ψ.ρ, dψ=dψ, d2ψ=d2ψ, k=k, δV=δV
        function integrand(r)
            δV(r,k)*ψ(r)*ψ(r)*(1/r^2)
        end
    end
    f
end

function ZeroOrderOperator(ψ1,ψ2)
    dψ = ∂r(ψ1);
    d2ψ = ∂r(dψ);
    f = let ψ=ψ1,ψ2=ψ2, s=ψ.s, l=ψ.l, ρ = ψ.ρ, dψ=dψ, d2ψ=d2ψ
        function integrand(r)
            Vr = (ρ^2*r^3)/(r-1) + l*(l+1) - (1-s^2)/r
            return ψ2(r)*(r*(r-1)*d2ψ(r) + dψ(r) - Vr*ψ(r))/r^2
        end
    end
    f
end


function FrequencyDerivativeOfOperator(ψs)
    f = let ψ=ψs, ρ = ψ.ρ
        function integrand(r)
            return (2*(im*ρ)*r^3*ψ(r))/(r^2*(r-1))
        end
    end
    f
end
