def design_liquid_crystal():
    """
    This script outlines the design principles for a single benzene ring-based
    liquid crystal material with transition temperatures around room temperature.
    """
    
    # Starting parameter for the alkyl chain length, as per suggestion F.
    n = 5
    
    print("--- Design Brief for a Room-Temperature Liquid Crystal ---")
    
    # A. Key requirements
    print("\n[A] Key Requirements:")
    print("- Core Structure: Single benzene ring")
    print("- Target Phases: Nematic or Smectic")
    print("- Target Transition Temperatures: Near room temperature (~20-25°C)")
    
    # B, C, D, E. Suggested Molecular Design
    print("\n[B, C, D, E] Suggested Molecular Design & Structure:")
    print("- A single benzene ring (Ph) is used as the rigid core.")
    print("- Functional groups are attached at the para positions (1 and 4).")
    print("- Feature 1: A flexible alkyl chain (C_n H_{2n+1}-) for phase behavior.")
    print("- Feature 2: A polar group (-CN) to enhance liquid crystal properties.")
    print(f"- Specific Example (for n={n}): 4-pentylbenzonitrile")
    
    # E. The general structure "equation" and its numbers
    print("\n[E] General Structure Equation:")
    print("The general chemical structure is represented by the formula:")
    print("C_n H_{2n+1} - Ph - CN")
    print("\nAnalyzing the numbers within this equation for the alkyl group (C_n H_{2n+1}):")
    
    # Explicitly printing the numbers from the equation as requested
    coefficient_of_n = 2
    constant_term = 1
    
    print(f"1. The number '{coefficient_of_n}' is the coefficient for 'n' in the hydrogen count.")
    print(f"2. The number '{constant_term}' is the constant added to calculate the hydrogen count.")
    
    hydrogen_count = coefficient_of_n * n + constant_term
    print(f"\nFor our starting example where n = {n}:")
    print(f"  - The chemical formula for the alkyl group is C{n}H({coefficient_of_n}*{n}+{constant_term}) = C{n}H{hydrogen_count}")

    # F. Tuning Strategy
    print("\n[F] Strategy for Tuning Transition Temperature:")
    print(f"- Start with a chain of length n={n} (pentyl).")
    print("- If the material's transition temperature is too high, you should INCREASE the chain length (n).")
    print("- If the transition temperature is too low, you should DECREASE the chain length (n).")
    print("- For fine-tuning, lateral substituents (e.g., -F, -CH3) can be added to the benzene ring.")

# Execute the design function
design_liquid_crystal()