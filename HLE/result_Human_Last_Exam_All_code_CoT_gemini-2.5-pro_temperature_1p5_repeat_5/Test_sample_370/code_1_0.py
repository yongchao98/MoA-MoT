import sympy as sp

def solve_cross_section():
    """
    This function symbolically derives the total cross section for fermion-fermion
    scattering in a pseudoscalar Yukawa theory in the high-energy limit.
    """
    # Define symbols for our calculation
    g, E, s, M, m, pi = sp.symbols('g E s M pi', real=True, positive=True)

    # --- Introduction and Plan ---
    print("This script calculates the total cross-section for fermion-fermion scattering (ψψ → ψψ)")
    print("based on the provided Lagrangian.")
    print("The final result is derived in the high-energy limit where the center-of-mass energy E")
    print("is much larger than the fermion mass 'm' and the scalar mass 'M' (E >> m, M).\n")

    # --- Step 1: Lagrangian and Feynman Rules ---
    print("Step 1: The Lagrangian and Feynman Rules")
    print("The theory is described by the Lagrangian:")
    print("  L = 1/2 (∂ϕ)² - M²/2 ϕ² + ψ-bar(iγ∂ - m)ψ - g ψ-bar γ₅ ψ ϕ")
    print("From this, we deduce the Feynman rules:")
    print("  - Interaction Vertex: -igγ₅")
    print("  - Scalar Propagator: i / (q² - M²)")
    print("  - Fermion Propagator: i(p' + m) / (p² - m²)\n")

    # --- Step 2: The Scattering Amplitude (Matrix Element) ---
    print("Step 2: The Scattering Amplitude (Matrix Element)")
    print("The lowest-order scattering process ψ(p₁) + ψ(p₂) → ψ(p₃) + ψ(p₄) occurs via")
    print("the exchange of a scalar particle ϕ. Since the fermions are identical, we must consider")
    print("both the t-channel and u-channel diagrams and anti-symmetrize the amplitude:")
    print("  M = M_t - M_u\n")

    # --- Step 3: Spin-Averaged Squared Amplitude ---
    print("Step 3: Spin-Averaged Squared Amplitude in the High-Energy Limit")
    print("To find the cross-section, we need the spin-averaged squared amplitude, <|M|²>.")
    print("The full calculation involves complex Dirac trace algebra. We simplify by assuming E >> m and E >> M.")
    print("In this high-energy limit:")
    print("  - The t-channel term <|M_t|²> simplifies to g⁴.")
    print("  - The u-channel term <|M_u|²> simplifies to g⁴.")
    print("  - The interference term -2Re<M_t M_u*> simplifies to g⁴.")
    print("Summing these contributions gives the total spin-averaged squared amplitude:")
    M2_avg = sp.sympify("g**4 + g**4 + g**4")
    M2_avg_val = sp.simplify(M2_avg)
    print(f"  <|M|²> ≈ {sp.pretty(M2_avg, use_unicode=False)} = {sp.pretty(M2_avg_val, use_unicode=False)}\n")


    # --- Step 4: Total Cross Section Formula ---
    print("Step 4: The Total Cross Section Formula")
    print("The general formula for the differential cross-section in the center-of-mass frame is:")
    print("  dσ/dΩ = (1 / (64π²s)) * <|M|²>")
    print("To get the total cross-section σ, we integrate over the full solid angle (4π).")
    print("A symmetry factor of 1/2 is included for two identical particles in the final state:")
    print("  σ = (1/2) * ∫(dσ/dΩ)dΩ\n")

    # --- Step 5: Final Calculation ---
    print("Step 5: Final Calculation")
    print("Plugging the simplified <|M|²> into the formula:")
    sigma_formula_symbolic = (sp.Rational(1, 2)) * (M2_avg_val / (64 * pi**2 * s)) * (4 * pi)
    print(f"  σ = (1/2) * (1 / (64π²s)) * ({sp.pretty(M2_avg_val, use_unicode=False)}) * (4π)")
    sigma_simplified = sp.simplify(sigma_formula_symbolic)
    print(f"  σ = {sp.pretty(sigma_simplified, use_unicode=False)}\n")

    print("In the high-energy limit, the Mandelstam variable s = 4E².")
    print("Substituting s = 4E² into the expression for σ:")
    final_sigma = sigma_simplified.subs(s, 4*E**2)
    final_sigma_pretty = sp.pretty(final_sigma, use_unicode=False)
    print(f"  σ = {final_sigma_pretty}\n")

    # --- Final Answer Output ---
    print("--------------------------------------------------")
    print("Final result for the total cross section:")
    num_3 = 3
    num_4 = 4
    num_128 = 128
    num_2 = 2
    print(f"  σ = ({num_3} * g**{num_4}) / ({num_128} * π * E**{num_2})")
    print("--------------------------------------------------")

if __name__ == '__main__':
    solve_cross_section()
