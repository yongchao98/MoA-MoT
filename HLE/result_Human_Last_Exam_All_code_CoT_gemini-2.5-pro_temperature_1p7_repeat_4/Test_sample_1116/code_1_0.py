def solve_chemistry_problem():
    """
    Analyzes the reaction and NMR data to identify the unknown product.
    Prints the explanation and the structure of the identified compound.
    """
    # Define molecular fragments for clarity
    # R represents the 2-ethylhexyl group: -CH2CH(CH2CH3)(CH2)3CH3
    # CORE represents the central dithienoisoindole-dione unit.
    R = "CH2CH(Et)Bu"
    CORE = "[DTI-Core]"
    
    # Represent the pendant thiophene arms
    unreacted_arm = f"({R})-Thiophene-"
    reacted_arm = f"({R})-Thiophene(Br)-"
    
    # Define the molecules in the reaction
    starting_material = f"{unreacted_arm}{CORE}{unreacted_arm}"
    mono_bromo_product = f"{reacted_arm}{CORE}{unreacted_arm}"
    
    print("Based on the reaction and NMR data, the new spot is the mono-brominated product.\n")
    print("### Reaction Analysis ###\n")
    print("The reaction that occurred is an incomplete bromination:\n")
    print(f"Starting Material: {starting_material}")
    print(f"Product Identified: {mono_bromo_product}\n")
    
    print("### Structural Explanation ###\n")
    print("The starting material is symmetrical. The intended product was the di-brominated compound, which would also be symmetrical.")
    print("However, the isolated product is the result of bromination on only ONE of the two thiophene rings, making the entire molecule asymmetrical.\n")

    print("### H-NMR Data Correlation ###\n")
    print("This asymmetry explains why 'three peaks that are larger than 6.0 ppm' are observed in the H-NMR spectrum.")
    print("These three peaks correspond to the three chemically non-equivalent protons on the two pendant thiophene rings:\n")

    protons = {
        "Peak 1": "The proton at position 5 of the UN-REACTED thiophene ring.",
        "Peak 2": "The proton at position 3 of the UN-REACTED thiophene ring.",
        "Peak 3": "The proton at position 3 of the REACTED (brominated) thiophene ring."
    }
    
    # The prompt requested to "output each number in the final equation".
    # We interpret this as identifying the 3 protons corresponding to the 3 NMR peaks.
    print("Final Equation (Proton origins):")
    for peak, description in protons.items():
        print(f" - {peak}: {description}")
        
    print("\nVisual representation of the protons on the two different thiophene arms of the product:")
    print("""
    Arm 1 (Reacted):                      Arm 2 (Unreacted):
    
          Br                                     H (Peak 1)
          |                                      |
    //===\\ /                                  //===\\ /
 R--C    C -- CORE -- ...               ...-- CORE -- C    C -- R
    \\ S //                                         \\ S //
     |                                               |
     H (Peak 3)                                    H (Peak 2)
    """)

# Execute the function to print the solution
solve_chemistry_problem()
