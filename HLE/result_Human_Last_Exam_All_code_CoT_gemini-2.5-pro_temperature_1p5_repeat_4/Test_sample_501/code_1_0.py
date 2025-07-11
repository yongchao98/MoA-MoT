import math

def display_force_law_derivation():
    """
    This function uses print statements to display the derivation and the final
    force law for a thermally isolated freely jointed polymer chain.
    """
    print("Derivation of the Force Law for a Thermally Isolated Polymer Chain:")
    print("="*60)

    # Step 1: State the physical principle
    print("1. For a thermally isolated (adiabatic) and reversible process, the total entropy is constant.")
    print("   S_total = S_kinetic(E_k) + S_configurational(x) = Constant\n")

    # Step 2: Define the entropy components
    print("2. The entropy components for a chain with n segments of length l are:")
    print("   S_kinetic(E_k) = (3n/2)*k_B*ln(E_k) + C1")
    print("   S_config(x)  = C2 - (3*k_B*x^2) / (2*n*l^2)\n")

    # Step 3: Find the energy-extension relation
    print("3. By setting S_total(E_k(x), x) = S_total(E(0), 0), we find how energy depends on extension:")
    print("   (3n/2)*k_B*ln(E_k(x)) - (3*k_B*x^2)/(2*n*l^2) = (3n/2)*k_B*ln(E(0))")
    print("   Solving for E_k(x) gives:")
    print("   E_k(x) = E(0) * exp(x^2 / (n^2 * l^2))\n")

    # Step 4: Derive the force equation
    print("4. The force F is found using F = T*(dS_config/dx) and T = 2*E_k/(3*n*k_B):")
    print("   F(x) = (2*E_k(x)/(3*n*k_B)) * ( - (3*k_B*x) / (n*l^2) )")
    print("   F(x) = - (2 * E_k(x) * x) / (n^2 * l^2)\n")

    # Step 5: Final Result
    print("5. Substituting E_k(x) from step 3 into the force equation from step 4 yields the final law.")
    print("="*60)
    print("Final Force Law:")
    print("F(x) = - (2 * E(0) / (n^2 * l^2)) * x * exp(x^2 / (n^2 * l^2))")
    print("="*60)
    
    print("As requested, here are the numerical constants and exponents in the final equation:")
    print("\nThe equation has the form: F(x) = - ( Pre-factor ) * x * exp( Exponent_Term )")
    
    print("\nIn the Pre-factor '(2 * E(0) / (n^2 * l^2))':")
    print(f"  - The numerical multiplier is: 2")
    print(f"  - The exponent on 'n' is: 2")
    print(f"  - The exponent on 'l' is: 2")

    print("\nIn the main term ' * x':")
    print(f"  - The exponent on 'x' is: 1")

    print("\nIn the Exponent Term 'exp(x^2 / (n^2 * l^2))':")
    print(f"  - The exponent on 'x' inside the exp() is: 2")
    print(f"  - The exponent on 'n' inside the exp() is: 2")
    print(f"  - The exponent on 'l' inside the exp() is: 2")


# Execute the function to print the explanation and result.
display_force_law_derivation()