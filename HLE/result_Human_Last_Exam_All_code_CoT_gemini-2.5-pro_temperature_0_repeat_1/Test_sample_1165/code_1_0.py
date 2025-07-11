import math

def solve_ode_fluctuations():
    """
    This function provides a step-by-step derivation for the fluctuation magnitude R(epsilon)
    and discusses the scaling for a different distribution of z_i.
    """

    print("### Part 1: Estimating R(epsilon) for Uniformly Distributed z_i ###")
    print("\nThe ODE is y'' - ε*y' = ε² * Σδ(x - zᵢ), with y(0)=1, y(ε⁻¹)=0.")
    
    print("\nStep 1: Asymptotic Expansion of the ODE")
    print("We expand y(x) = y₀(x) + ε*y₁(x) + ε²*y₂(x) + ...")
    print("Substituting this into the ODE and collecting terms by powers of ε gives:")
    print("O(1): y₀'' = 0")
    print("O(ε): y₁'' - y₀' = 0")
    print("O(ε²): y₂'' - y₁' = Σδ(x - zᵢ)")
    print("\nApplying boundary conditions y(0)=1, y(ε⁻¹)=0 to each order (yᵢ(0)=yᵢ(ε⁻¹)=0 for i>0), we solve for the first two terms:")
    print("y₀(x) = 1 - εx")
    print("y₁(x) = (x/2) * (1 - εx)")
    print("These terms are deterministic.")

    print("\nStep 2: Relating R to the Variance of y₂(x)")
    print("The fluctuation relative to the leading order term is y(x) - y₀(x) ≈ ε*y₁(x) + ε²*y₂(x).")
    print("The randomness comes from y₂(x), as y₀ and y₁ are deterministic.")
    print("The definition of R is R = (max_x Var[y(x) - y₀(x)])^(1/2).")
    print("R² = max_x Var[y(x) - y₀(x)] ≈ max_x Var[ε²*y₂(x)] = ε⁴ * max_x Var[y₂(x)]")
    
    print("\nStep 3: Solving for the Stochastic Part of y₂(x)")
    print("The equation for y₂ is y₂'' = y₁' + Σδ(x - zᵢ).")
    print("The variance of y₂ only depends on its stochastic part, y₂_stoch, where y₂_stoch'' = Σδ(x - zᵢ).")
    print("Using the Green's function G(x, ξ) for the operator d²/dx² with zero boundary conditions, we get:")
    print("y₂_stoch(x) = Σ G(x, zᵢ)")
    print("where G(x, ξ) = x(1 - εξ) for ξ > x and ξ(1 - εx) for ξ < x.")

    print("\nStep 4: Calculating Var[y₂(x)]")
    print("Var[y₂(x)] = Var[Σ G(x, zᵢ)].")
    print("A key property of order statistics (zᵢ) is that the distribution of Σf(zᵢ) is the same as Σf(Uᵢ) where Uᵢ are i.i.d.")
    print("So, Var[Σ G(x, zᵢ)] = Var[Σ G(x, Uᵢ)] where Uᵢ are i.i.d. Uniform on [0, ε⁻¹].")
    print("Since Uᵢ are independent, Var[Σ G(x, Uᵢ)] = N * Var[G(x, U)], where N = ε⁻¹ - 1 ≈ ε⁻¹.")
    print("We calculate Var[G(x, U)] = E[G(x, U)²] - (E[G(x, U)])².")
    print("After integration, we find Var[G(x, U)] = (1/12) * x² * (1 - εx)².")
    print("Therefore, Var[y₂(x)] ≈ ε⁻¹ * (1/12) * x² * (1 - εx)².")

    print("\nStep 5: Finding the Maximum Fluctuation R(ε)")
    print("R² ≈ ε⁴ * max_x [ε⁻¹ * (1/12) * x² * (1 - εx)²]")
    print("R² ≈ (ε³/12) * max_x [x² * (1 - εx)²]")
    print("To maximize, let u = εx. The term to maximize is (u/ε)²(1-u)² = ε⁻²u²(1-u)². The domain for u is [0, 1].")
    print("R² ≈ (ε³/12) * ε⁻² * max_u [u²(1-u)²]")
    print("R² ≈ (ε/12) * max_u [u²(1-u)²]")
    print("The maximum of u(1-u) is 1/4 at u=1/2. So the max of u²(1-u)² is (1/4)² = 1/16.")
    c_squared_inv = 12.0 * 16.0
    c = 1.0 / math.sqrt(c_squared_inv)
    sqrt3 = math.sqrt(3)
    print(f"R² ≈ (ε/12) * (1/16) = ε / 192.")
    print("The final estimate for R(ε) is:")
    print(f"R(ε) = (1 / sqrt({c_squared_inv:.0f})) * ε^(1/2) = (1 / (8 * sqrt(3))) * ε^(1/2)")
    print(f"R(ε) ≈ {c:.4f} * ε^(0.5)")

    print("\n\n### Part 2: Scaling for Normally Distributed z_i ###")
    print("If zᵢ ~ Normal(i, 0.5) are i.i.d., the situation changes.")
    print("The zᵢ are independent, so Var[y₂(x)] = Σ Var[G(x, zᵢ)].")
    print("The variance of zᵢ (σ²=0.25) is small. We can approximate Var[G(x, zᵢ)] ≈ (∂G/∂z(x, i))² * σ².")
    print("The derivative ∂G/∂z is O(1). So Var[G(x, zᵢ)] is O(1).")
    print("Approximating the sum over N ≈ ε⁻¹ terms by an integral: Var[y₂(x)] ≈ σ² * ∫[0, ε⁻¹] (∂G/∂z(x, u))² du.")
    print("The integral evaluates to (x/ε)(1 - εx), which is O(ε⁻¹).")
    print("So, Var[y₂(x)] ≈ σ² * O(ε⁻¹).")
    print("Then R² = ε⁴ * max_x Var[y₂(x)] ≈ ε⁴ * O(ε⁻¹) = O(ε³).")
    print("This gives R(ε) ~ ε^(3/2).")
    print("\nConclusion: The scaling for R(ε) is expected to change from ε^(1/2) to ε^(3/2).")

if __name__ == '__main__':
    solve_ode_fluctuations()