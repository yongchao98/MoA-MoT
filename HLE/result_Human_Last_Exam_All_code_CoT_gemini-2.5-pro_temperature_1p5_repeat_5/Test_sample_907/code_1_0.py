import sympy as sp

def display_absorption_equations():
    """
    This function uses the sympy library to symbolically construct and print
    the equations for the absorption cross-section in two different cases.
    """

    # --- Define all symbolic variables used in the equations ---
    C = sp.Symbol('C', positive=True, real=True) # Proportionality constant
    omega_L = sp.Symbol('ω_L', real=True)         # Laser carrier frequency
    omega_0 = sp.Symbol('ω_0', positive=True, real=True) # Molecular transition frequency
    mu_e0 = sp.Symbol('μ_e0', real=True)          # Molecular transition dipole moment
    tau = sp.Symbol('τ', positive=True, real=True)         # Pulse duration
    N = sp.Symbol('N', integer=True, positive=True) # Number of molecules
    J = sp.Symbol('J', real=True)                  # Near-neighbor coupling energy
    hbar = sp.Symbol('ħ', positive=True, real=True) # Reduced Planck constant
    j = sp.Symbol('j', integer=True)               # Exciton state index
    sigma_a = sp.Symbol('σ_a')
    sigma_b = sp.Symbol('σ_b')

    # --- Print Header and Symbol Definitions ---
    print("--------------------------------------------------------------------")
    print("Equations for Absorption Cross-Section")
    print("--------------------------------------------------------------------\n")
    print("This script provides the equations for the absorption cross-section σ(ω_L)")
    print("for a chain of molecules absorbing an ultrashort Gaussian-shape laser pulse.\n")
    print("Symbols used:")
    sp.pprint(sp.Eq(C, sp.Symbol("Proportionality constant")), use_unicode=True)
    sp.pprint(sp.Eq(omega_L, sp.Symbol("Laser carrier frequency")), use_unicode=True)
    sp.pprint(sp.Eq(omega_0, sp.Symbol("Single molecule transition frequency (E_exc - E_gnd)/ħ")), use_unicode=True)
    sp.pprint(sp.Eq(mu_e0, sp.Symbol("Single molecule transition dipole moment")), use_unicode=True)
    sp.pprint(sp.Eq(tau, sp.Symbol("Gaussian pulse duration")), use_unicode=True)
    sp.pprint(sp.Eq(N, sp.Symbol("Number of molecules in the chain")), use_unicode=True)
    sp.pprint(sp.Eq(J, sp.Symbol("Near-neighbor interaction energy")), use_unicode=True)
    sp.pprint(sp.Eq(hbar, sp.Symbol("Reduced Planck constant")), use_unicode=True)
    sp.pprint(sp.Eq(j, sp.Symbol("Index for exciton states (j=1, 2, ...)")), use_unicode=True)
    print("\n" + "="*70 + "\n")

    # --- Case a) No interaction between molecules ---
    print("Case a) The interaction between molecules can be neglected.\n")
    print("In this scenario, molecules are independent. The total absorption is the sum of")
    print("the absorptions of N individual molecules. This creates a single absorption")
    print("peak at the molecular transition frequency ω_0. The peak's Gaussian shape")
    print("is determined by the laser pulse's spectral profile.\n")
    print("The equation for the total absorption cross-section is:")
    
    eq_a = C * N * omega_0 * mu_e0**2 * sp.exp(-(omega_0 - omega_L)**2 * tau**2 / 2)
    sp.pprint(sp.Eq(sigma_a(omega_L), eq_a), use_unicode=True)
    print("\n" + "="*70 + "\n")

    # --- Case b) Near-neighbor interaction is considered ---
    print("Case b) The interaction between near-neighbors should be considered.\n")
    print("Here, interactions lead to delocalized exciton states forming an energy band.")
    print("For a linear chain model, only specific exciton states (those with odd index j)")
    print("are optically active. The spectrum is a series of peaks, each corresponding")
    print("to a transition to one of these allowed exciton states.\n")

    # Energy of the j-th exciton state relative to the ground state
    omega_j_expr = omega_0 + (2 * J / hbar) * sp.cos(sp.pi * j / (N + 1))
    print("The transition frequency for the j-th exciton state is:")
    sp.pprint(sp.Eq(sp.Symbol('ω_j'), omega_j_expr), use_unicode=True)
    print("")

    # The square of the transition dipole moment to the j-th exciton state (for odd j)
    mu_j_sq_expr = mu_e0**2 * (2 / (N + 1)) * (sp.cot(sp.pi * j / (2 * (N + 1))))**2
    print("The squared transition dipole moment for an optically active state (odd j) is:")
    sp.pprint(sp.Eq(sp.Symbol('|μ_j|')**2, mu_j_sq_expr), use_unicode=True)
    print("")

    # The total absorption cross-section is the sum over all optically active modes (j = 1, 3, 5, ...)
    print("The total absorption cross-section is the sum over all optically active modes (j=1, 3, 5, ...):")
    final_expression_b = (
        f"σ_b(ω_L) = C * Σ_{{j=1,3,5,...}}^{{N}} [ "
        f"ω_j * |μ_j|² * exp(-(ω_j - ω_L)² * τ² / 2) ]"
    )
    print(final_expression_b)
    print("\nSubstituting the expressions for ω_j and |μ_j|² gives the full formula.")
    print("The strongest absorption typically corresponds to the j=1 exciton state.")

if __name__ == '__main__':
    display_absorption_equations()