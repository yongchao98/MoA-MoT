def solve_sphere_energy_dissipation():
    """
    This function provides a step-by-step derivation for the Joule heat dissipated
    by a shrinking, charge-leaking sphere and prints the final formula.
    """

    # Define symbolic representations for the physical quantities
    a = 'a'
    V = 'V'
    pi = 'π'
    epsilon_0 = 'ε₀'

    print("Derivation for Joule Heat of a Shrinking, Leaking Sphere:\n")

    print("Step 1: Initial Electrostatic Energy")
    print(f"A sphere with radius r={a} and potential V={V} has an initial charge Q_i = 4⋅{pi}⋅{epsilon_0}⋅{a}⋅{V}.")
    print(f"The initial stored electrostatic energy is U_i = (1/2)⋅Q_i⋅V = 2⋅{pi}⋅{epsilon_0}⋅{a}⋅{V}².")
    print("-" * 60)

    print("Step 2: Energy Conservation Principle")
    print("The initial energy U_i is converted into Joule heat (Q_Joule) and mechanical work (W_field).")
    print("The exact distribution depends on the process path (the relation between charge and radius).")
    print("-" * 60)
    
    print("Step 3: Assumption of Constant Surface Charge Density")
    print("To find a unique solution, we assume a physically plausible process where the surface charge density 'σ' remains constant.")
    print(f"This implies that at any radius r, the potential is V(r) = V⋅(r/a).")
    print("The process ends when the charge is zero, which corresponds to a final radius of r=0.")
    print("-" * 60)

    print("Step 4: Integration of Dissipated Heat")
    print("The incremental Joule heat is dQ_Joule = -V(r)⋅dq.")
    print("By integrating dQ_Joule from the start (r=a) to the end (r=0) under the constant 'σ' assumption, we get the total heat.")
    print("The integral evaluates to: Q_Joule = ∫[a,0] -(8⋅π⋅ε₀⋅V²/a²)⋅r² dr")
    print("-" * 60)

    print("Step 5: The Final Result")
    print("After performing the integration, the total Joule heat dissipated into the atmosphere is:")
    
    # Define the numerical coefficients to fulfill the user request
    numerator = 8
    denominator = 3
    
    # Print the final equation with clear formatting
    print(f"\n    Q_Joule = ({numerator}/{denominator}) ⋅ {pi} ⋅ {epsilon_0} ⋅ {a} ⋅ {V}²\n")

solve_sphere_energy_dissipation()