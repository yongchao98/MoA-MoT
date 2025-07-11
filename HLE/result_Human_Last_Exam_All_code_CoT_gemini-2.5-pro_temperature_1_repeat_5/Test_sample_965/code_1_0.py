import sympy

def solve_photon_rate():
    """
    Calculates and explains the rate of photon production for a two-level atom
    in a cavity, based on the Purcell effect.
    """

    # Define symbolic variables for clarity in the output
    g = sympy.Symbol('g')
    gamma_c = sympy.Symbol('γ_c')
    kappa = sympy.Symbol('κ')

    print("Step 1: The rate of photon production is given by the Purcell-enhanced spontaneous emission rate, Γ.")
    print("The standard formula in the weak coupling limit is Γ = 2 * g² / κ")
    print("where g is the atom-cavity coupling frequency and κ is the cavity energy decay rate.")
    print("-" * 20)

    print("Step 2: We interpret the variable γ_c from the problem's answer choices.")
    print("It's common for γ_c to denote the cavity field *amplitude* decay rate, which is half the energy decay rate.")
    print("Therefore, we assume the relation: κ = 2 * γ_c.")
    print("-" * 20)
    
    print("Step 3: Substitute κ = 2 * γ_c into the standard formula.")
    # Standard formula
    rate_formula = 2 * g**2 / kappa
    # Substitute kappa
    final_rate = rate_formula.subs(kappa, 2 * gamma_c)

    print(f"Γ = 2 * g² / κ  =>  Γ = 2 * g² / (2 * γ_c)  =>  Γ = {final_rate}")
    print("-" * 20)

    print("Step 4: Output the components of the final equation as requested.")
    # The final equation is g**2 / gamma_c
    numerator_coeff = 1
    numerator_vars = "g^2"
    denominator_coeff = 1
    denominator_vars = "γ_c"
    
    print(f"Numerator: {numerator_coeff} * {numerator_vars}")
    print(f"Denominator: {denominator_coeff} * {denominator_vars}")
    print(f"Final Rate Γ = ({numerator_vars}) / ({denominator_vars})")

solve_photon_rate()