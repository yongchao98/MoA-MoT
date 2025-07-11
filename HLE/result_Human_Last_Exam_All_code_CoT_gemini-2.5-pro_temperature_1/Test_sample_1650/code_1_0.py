import math

def derive_overlap_integral():
    """
    Prints the step-by-step derivation of the overlap integral for two 2s orbitals.
    """
    print("Derivation of the Overlap Integral S(2s, 2s) for H₂⁺\n")
    print("="*60)

    # Step 1: Definition of the integral
    print("\nStep 1: The Overlap Integral Definition")
    print("The overlap integral S between two atomic orbitals ψ_A and ψ_B, centered on nuclei A and B, is:")
    print("S = ∫ ψ_2s(r_A) * ψ_2s(r_B) dτ")
    print("where r_A and r_B are the distances of the electron from nucleus A and B, respectively.\n")

    # Step 2: The 2s wavefunction
    print("Step 2: The Hydrogen-like 2s Wavefunction")
    print("The normalized hydrogen-like 2s wavefunction (in atomic units) is given by:")
    print("ψ_2s(r) = N * (2 - ζr) * e^(-ζr/2)")
    print("where ζ is the effective nuclear charge and N = (ζ³/ (32π))^(1/2) is the normalization constant.\n")

    # Step 3: Transformation to elliptical coordinates
    print("Step 3: Transformation to Elliptical Coordinates")
    print("For a two-center system, it's best to use elliptical coordinates (λ, μ, φ):")
    print("  λ = (r_A + r_B) / R,   (1 ≤ λ < ∞)")
    print("  μ = (r_A - r_B) / R,   (-1 ≤ μ ≤ 1)")
    print("  φ = azimuthal angle,  (0 ≤ φ ≤ 2π)")
    print("The volume element is dτ = (R³/8) * (λ² - μ²) dλ dμ dφ.")
    print("Let's define a dimensionless quantity ρ = ζR. The integral becomes:")
    print("S = (ρ³/128) ∫[1,∞] dλ ∫[-1,1] dμ ∫[0,2π] dφ * [4 - 2ρλ + (ρ²/4)(λ² - μ²)] * e^(-ρλ/2) * (λ² - μ²)\n")

    # Step 4: Integration over φ and μ
    print("Step 4: Integration over φ and μ")
    print("The integral over φ is trivial, yielding 2π.")
    print("The integral over μ, I_μ = ∫[-1,1] [...], results in a polynomial in λ:")
    print("I_μ(λ) = (ρ²/2)λ⁴ - 4ρλ³ + (8 - ρ²/3)λ² + (4ρ/3)λ + (ρ²/10 - 8/3)\n")

    # Step 5: Integration over λ
    print("Step 5: Integration over λ")
    print("The remaining integral over λ is:")
    print("S = (ρ³/128) ∫[1,∞] I_μ(λ) * e^(-ρλ/2) dλ")
    print("This is solved using auxiliary functions A_n(α) = ∫[1,∞] xⁿ * e^(-αx) dx, with α = ρ/2.")
    print("After substituting the solutions for A_n(α) and collecting terms, a significant simplification occurs.\n")

    # Step 6: Final result
    print("Step 6: Final Analytical Expression")
    print("After combining all parts, the final analytical expression for the overlap integral S is:")
    print("\n" + "-"*25 + " FINAL EQUATION " + "-"*25)
    
    # Printing the final equation with each number explicitly stated
    # The numbers in the equation are: 1, 1, 2, 1, 12, 1, 240
    print("\nThe overlap integral S as a function of ρ = ζ*R is:")
    print("S(ρ) = e^(-ρ / 2) * (1 + (1 / 2)*ρ + (1 / 12)*ρ^2 + (1 / 240)*ρ^4)")
    
    print("\nIn terms of the physical quantities ζ (effective nuclear charge) and R (internuclear distance):")
    print("S(ζ, R) = e^(-(ζ*R) / 2) * (1 + (1/2)*(ζ*R) + (1/12)*(ζ*R)^2 + (1/240)*(ζ*R)^4)")
    print("-" * 68)


if __name__ == '__main__':
    derive_overlap_integral()

<<<S(ζ, R) = e^(-(ζ*R) / 2) * (1 + (1/2)*(ζ*R) + (1/12)*(ζ*R)^2 + (1/240)*(ζ*R)^4)>>>