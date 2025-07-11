import sympy as sp

def analyze_soliton_stability():
    """
    Analyzes the stability of a 3D localized soliton based on the Hobart-Derrick scaling theorem.
    The Hamiltonian density is H = A*(grad m)^2 + D*m.(curl m).
    """

    # --- Introduction to the Scaling Argument ---
    print("Analyzing the stability of a localized soliton solution using a scaling argument (Hobart-Derrick theorem).")
    print("We assume a reference solution and scale its spatial dimensions by a factor lambda.")
    print("The total energy E of the scaled solution can then be expressed as a function of lambda.\n")
    
    # --- 1. Symbolic Energy Formulation ---
    print("--- 1. Symbolic Energy Formulation ---")
    l, E_A, E_D = sp.symbols('lambda E_A E_D')
    
    # In 3D, for a scaling r -> lambda*r, the energy terms scale as follows:
    # - The exchange energy integral, proportional to integral[(grad m)^2] dV, scales as (1/lambda^2) * lambda^3 = lambda.
    # - The DMI energy integral, proportional to integral[m.(curl m)] dV, scales as (1/lambda) * lambda^3 = lambda^2.
    
    E_total = E_A * l + E_D * l**2

    print(f"The total energy as a function of the soliton size lambda is: E(lambda) = E_A * lambda + E_D * lambda^2")
    print("Here, E_A and E_D are the exchange and DMI energy values for the unscaled (lambda=1) solution.\n")

    # --- 2. Conditions for a Stable Minimum ---
    print("--- 2. Conditions for a Stable Minimum ---")
    dE_dl = sp.diff(E_total, l)
    d2E_dl2 = sp.diff(dE_dl, l)

    print(f"To find a stationary point (potential minimum), we set the first derivative to zero: dE/dlambda = {dE_dl} = 0")
    print(f"For this point to be a stable minimum, the second derivative must be positive: d^2E/dlambda^2 = {d2E_dl2} > 0")

    print("\nEvaluating at the unscaled solution (lambda=1):")
    print(f"Condition 1 (stationary): {dE_dl.subs(l, 1)} = 0  =>  E_A = -2*E_D")
    print(f"Condition 2 (stability): {d2E_dl2} > 0             =>  E_D > 0")

    # --- 3. The Physical Contradiction ---
    print("\n--- 3. The Physical Contradiction ---")
    print("Let's analyze the physical constraints on E_A and E_D:")
    print("a) For a non-uniform soliton, its spatial gradients must be non-zero.")
    print("   The exchange energy E_A = integral[A*(grad m)^2] dV, with A > 0, is therefore strictly positive (E_A > 0).")
    print("b) The DMI drives the formation of the chiral soliton. For this to be energetically favorable, it must lower the total energy.")
    print("   This means the DMI energy E_D = integral[D*m.(curl m)] dV must be negative (E_D < 0).")

    print("\nNow we check for contradictions:")
    print(f"From physical constraint (a), we need E_A > 0.")
    print(f"From Condition 1 (E_A = -2*E_D), if E_A > 0, then -2*E_D must be positive, which means E_D must be negative (E_D < 0).")
    print(f"This aligns with physical constraint (b). So, a stationary point can exist if E_D is negative.")
    
    print(f"\nHowever, Condition 2 for stability requires E_D > 0.")
    print("This directly contradicts the requirement that E_D must be negative for a stationary point to exist.")
    print("\nConclusion: It is impossible to satisfy the conditions for both a stationary point and stability simultaneously.")
    print("Any stationary point is an energy maximum (saddle point), not a minimum.\n")

    # --- 4. Numerical Illustration ---
    print("--- 4. Numerical Illustration ---")
    # Choose values that satisfy the condition for an extremum: E_A > 0, E_D < 0.
    E_A_val = 4.0
    E_D_val = -1.0
    
    # Calculate the stationary point from dE/dl = E_A + 2*E_D*l = 0
    lambda_extremum = -E_A_val / (2 * E_D_val)

    print(f"Let's choose example values: E_A = {E_A_val}, E_D = {E_D_val}.")
    print("These values satisfy the physical requirement for a soliton to attempt to form (E_A > 0, E_D < 0).")
    print(f"The total energy equation is: E(lambda) = {E_A_val} * lambda + ({E_D_val}) * lambda^2")
    print(f"The energy has a stationary point at lambda = -({E_A_val}) / (2 * ({E_D_val})) = {lambda_extremum}")
    
    second_deriv_val = 2 * E_D_val
    print(f"The second derivative at this point is 2 * E_D = 2 * ({E_D_val}) = {second_deriv_val}")
    print("Since the second derivative is negative, the solution corresponds to an energy maximum and is UNSTABLE.")
    print("The soliton will either collapse (lambda -> 0) or expand without limit (lambda -> infinity) to lower its energy.")

# Execute the analysis
analyze_soliton_stability()