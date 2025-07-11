def display_force_law_derivation():
    """
    This function explains the derivation of the force law for a thermally
    isolated, freely jointed polymer chain and prints the final formula.
    """

    print("Derivation of the Force Law for a Thermally Isolated Polymer Chain:")
    print("="*70)
    print("Step 1: Define the System and Thermodynamic Principles")
    print("  - The system is a freely jointed chain with n segments of length l.")
    print("  - It is thermally isolated, so for a slow, reversible extension, the total entropy S is constant (dS = 0).")
    print("  - The total energy E (which is purely kinetic) changes only due to work done on it: dE = F dx.")
    print("\nStep 2: Decompose the Entropy")
    print("  - Total entropy S = S_conf(x) + S_kin(E).")
    print("  - Configurational entropy (for small extension x): S_conf(x) = C1 - (3 * k_B * x^2) / (2 * n * l^2).")
    print("  - Kinetic entropy (for f degrees of freedom): S_kin(E) = (f/2) * k_B * ln(E) + C2.")
    print("  - For a large chain with n masses, the number of degrees of freedom f is approximately 2n.")
    print("\nStep 3: Relate Force to Entropy")
    print("  - From dS = 0, we have (∂S/∂x)_E dx + (∂S/∂E)_x dE = 0.")
    print("  - This gives (dS_conf/dx) dx + (dS_kin/dE) dE = 0.")
    print("  - Substituting dE = F dx, we solve for the force F: F = - (dS_conf/dx) / (dS_kin/dE).")
    print("\nStep 4: Calculate the Force in Terms of Instantaneous Energy E(x)")
    print("  - dS_conf/dx = - (3 * k_B * x) / (n * l^2).")
    print("  - dS_kin/dE = (f/2) * k_B / E ≈ (n * k_B) / E.")
    print("  - F = - [ - (3 * k_B * x) / (n * l^2) ] / [ (n * k_B) / E ] = (3 * E * x) / (n^2 * l^2).")
    print("\nStep 5: Express Force in Terms of Initial Energy E(0)")
    print("  - We solve the differential equation dE/dx = F = (3 * E * x) / (n^2 * l^2).")
    print("  - Separating variables: dE/E = (3 * x / (n^2 * l^2)) dx.")
    print("  - Integrating from x=0 to x gives: E(x) = E(0) * exp( (3 * x^2) / (2 * n^2 * l^2) ).")
    print("  - Substituting this E(x) back into the force equation gives the final result.")
    print("="*70)
    print("\nFinal Force Law:")
    
    # The instruction "output each number in the final equation" is interpreted
    # as printing the symbolic formula clearly.
    
    equation = "F(x) = (3 * E(0) * x / (n^2 * l^2)) * exp( (3 * x^2) / (2 * n^2 * l^2) )"
    
    print(equation)
    
    print("\nWhere:")
    print("  F(x) is the force between the polymer ends.")
    print("  x    is the separation of the ends.")
    print("  E(0) is the kinetic energy of the polymer at zero extension (x=0).")
    print("  n    is the number of segments in the polymer chain (assumed to be large).")
    print("  l    is the length of each segment.")
    print("  exp  is the exponential function.")

# Execute the function to display the result.
display_force_law_derivation()