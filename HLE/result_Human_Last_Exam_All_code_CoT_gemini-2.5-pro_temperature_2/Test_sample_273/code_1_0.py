import sympy

def analyze_soliton_stability():
    """
    Analyzes soliton stability in a 3D system with exchange and DMI
    using the Hobart-Derrick scaling argument.
    """

    # --- Setup ---
    # Define symbolic variables for the energy terms (at lambda=1) and the scaling factor.
    # E_A: Exchange energy, must be positive for any non-uniform texture.
    # E_D: DMI energy.
    # lam: The dimensionless scaling factor lambda.
    E_A, E_D, lam = sympy.symbols('E_A E_D lambda')
    
    # Define spatial dimension
    d = 3

    print(f"--- Soliton Stability Analysis in {d}D ---")
    print("This script uses a scaling argument to test stability.")
    print("We scale a hypothetical soliton by a factor lambda and analyze the energy E(lambda).")
    
    # --- Energy Scaling ---
    # The scaling of each energy term depends on the power of the gradient and the volume element.
    # Exchange energy term scales as lambda^(d-2)
    # DMI energy term scales as lambda^(d-1)
    exchange_exponent = d - 2
    dmi_exponent = d - 1
    
    # Total energy as a function of the scaling parameter lambda
    total_energy = E_A * lam**exchange_exponent + E_D * lam**dmi_exponent
    
    print("\nThe total energy E(lambda) of the scaled soliton is:")
    print(f"E(lambda) = E_A * lambda**({d} - 2) + E_D * lambda**({d} - 1)")
    print(f"E(lambda) = E_A * lambda**{exchange_exponent} + E_D * lambda**{dmi_exponent}\n")

    # --- First Condition: Energy Extremum ---
    # For a stable solution at lambda=1, the first derivative of energy must be zero.
    dE_dlam = sympy.diff(total_energy, lam)
    
    print("--- Condition 1: Energy Extremum (dE/d(lambda) = 0 at lambda = 1) ---")
    print(f"The first derivative dE/d(lambda) is: {dE_dlam}")

    # Set lambda=1 to get the equilibrium condition (Virial Theorem)
    virial_theorem_eq = sympy.Eq(dE_dlam.subs(lam, 1), 0)
    print("The condition for an energy extremum at lambda=1 is (Virial Theorem):")
    # Pretty print the final equation with all numbers
    print(f"Equation: {virial_theorem_eq.lhs.args[1].args[0]}*E_A + {virial_theorem_eq.lhs.args[0].args[0]}*E_D = 0")


    # --- Second Condition: Stability ---
    # For the extremum to be a MINIMUM, the second derivative must be positive.
    d2E_dlam2 = sympy.diff(dE_dlam, lam)

    print("\n--- Condition 2: Stability (d^2E/d(lambda)^2 > 0 at lambda = 1) ---")
    print(f"The second derivative d^2E/d(lambda)^2 is: {d2E_dlam2}")
    
    second_derivative_at_1 = d2E_dlam2.subs(lam, 1)
    print("The stability condition at lambda=1 requires:")
    print(f"Equation: {second_derivative_at_1} > 0\n")

    # --- Analysis ---
    print("--- Analysis of the Conditions ---")
    print("From Condition 1 (Virial Theorem):")
    print(f"   {exchange_exponent}*E_A + {dmi_exponent}*E_D = 0")

    # To satisfy this, if E_A > 0 (cost of creating a texture), E_D must be negative.
    # We can express E_D in terms of E_A from this condition.
    e_d_solution = sympy.solve(virial_theorem_eq, E_D)[0]
    print(f"   For a non-trivial soliton (E_A > 0), this means E_D = {e_d_solution}")

    # Now substitute this into the second condition (stability)
    stability_check = second_derivative_at_1.subs(E_D, e_d_solution)

    print("\nNow we check Condition 2 (Stability) using the result from Condition 1:")
    print(f"   Value of second derivative = {second_derivative_at_1}")
    print(f"   Substituting E_D = {e_d_solution}, we get: {stability_check}")

    print("\nSince E_A (exchange energy of a texture) must be positive, the second derivative is always negative.")
    print("A negative second derivative indicates that the energy is at a MAXIMUM, not a minimum.")
    
    # --- Final Conclusion ---
    print("\n--- Conclusion ---")
    print("Any localized non-trivial configuration is unstable. It will tend to either collapse to a point (lambda -> 0) or expand indefinitely (lambda -> infinity).")
    print("Therefore, it is not possible to stabilize a localized soliton with only these two energy terms in 3D.")


if __name__ == '__main__':
    analyze_soliton_stability()
    print("\n\n<<<No>>>")