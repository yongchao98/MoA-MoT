import sympy as sp

def analyze_soliton_stability():
    """
    Performs a scaling analysis (Hobart-Derrick theorem) on a 3D magnetic
    Hamiltonian with exchange and DMI terms to determine soliton stability.
    """
    # Define symbols for the analysis
    # l (lambda) is the dimensionless scaling factor for the soliton's size.
    # c1 represents the integrated exchange energy of the initial soliton (must be positive).
    # c2 represents the integrated DMI energy of the initial soliton (must be negative to favor a soliton).
    l, c1, c2 = sp.symbols('lambda c1 c2', real=True, positive=True)
    
    # In our physical model, the DMI term must lower the energy for a chiral soliton
    # to form, so its coefficient is negative. Let's represent the DMI contribution
    # as -c2 where c2 is a positive value.
    E_DMI = -c2 * l**2
    E_exchange = c1 * l

    # The total energy E is the sum of the exchange and DMI contributions.
    E = E_exchange + E_DMI

    # --- Step 1: Explain the model and energy scaling ---
    print("--- Analysis of Soliton Stability in 3D ---")
    print("Hamiltonian Density: H = A*(grad m)^2 + D*m.(curl m)")
    print("\nWe use a scaling argument. Let a soliton have a characteristic size L.")
    print("We scale the size by a factor lambda: L -> lambda * L.")
    print("\nScaling of Energy Contributions:")
    print("1. Heisenberg Exchange Energy (resists gradients): E_exchange ~ A * L")
    print("   In our model: E_exchange(lambda) = c1 * lambda, with c1 > 0.")
    print("2. Dzyaloshinskii-Moriya Energy (favors chirality): E_DMI ~ D * L^2")
    print("   In our model: E_DMI(lambda) = -c2 * lambda^2, with c2 > 0.")
    
    # --- Step 2: Formulate and analyze the total energy function ---
    print("\nThe total energy as a function of the scaling factor lambda is:")
    print(f"E(lambda) = {c1}*lambda - {c2}*lambda**2")

    # --- Step 3: Find the energy extremum ---
    print("\nTo find a size where the energy is at an extremum (min or max),")
    print("we take the first derivative of E(lambda) and set it to zero.")
    
    dE_dl = sp.diff(E, l)
    print(f"First derivative: dE/d(lambda) = {dE_dl}")
    
    # Solve for the lambda value at the extremum
    extremum_eq = sp.Eq(dE_dl, 0)
    l_extremum_sol = sp.solve(extremum_eq, l)
    l_extremum = l_extremum_sol[0]
    
    print("\nSetting the derivative to zero gives the final equation for the extremum:")
    print(f"{dE_dl} = 0")
    print("The numbers/coefficients in this equation are: 1 (for c1), -2 (for c2), and 0.")
    print(f"\nThis gives a potential equilibrium size at: lambda_extremum = {l_extremum}")
    print("Since c1 and c2 are positive, this corresponds to a real, finite-sized solution.")

    # --- Step 4: Check the stability of the extremum ---
    print("\nTo determine if this is a stable minimum, we check the second derivative.")
    d2E_dl2 = sp.diff(dE_dl, l)
    print(f"Second derivative: d^2E/d(lambda)^2 = {d2E_dl2}")
    
    print("\nSince c2 is a positive constant, the second derivative is always negative.")
    print("A negative second derivative means the energy is at a MAXIMUM at this size.")

    # --- Step 5: Final Conclusion ---
    print("\n--- Conclusion ---")
    print("The analysis shows that any localized soliton solution corresponds to an")
    print("energy maximum, not a minimum. The soliton is therefore unstable.")
    print("It can lower its energy by either collapsing to zero size (lambda -> 0)")
    print("or by expanding indefinitely (lambda -> infinity).")
    print("\nTherefore, it is not possible to stabilize a localized soliton with only")
    print("the given exchange and DMI terms in 3D.")

analyze_soliton_stability()