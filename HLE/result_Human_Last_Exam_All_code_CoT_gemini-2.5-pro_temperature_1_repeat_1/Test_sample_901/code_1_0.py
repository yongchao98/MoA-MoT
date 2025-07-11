def identify_elimination_product():
    """
    Analyzes the reaction of (1S,2R)-1-bromo-2-methylcyclohexane with potassium tert-butoxide
    and identifies the major product.
    """
    
    # Define the reaction components
    substrate = "(1S,2R)-1-bromo-2-methylcyclohexane"
    reagent = "potassium tert-butoxide"
    
    # Print the analysis steps
    print("### Reaction Analysis ###\n")
    print(f"Substrate: {substrate}")
    print(f"Reagent: {reagent}\n")
    
    print("--- Step 1: Determine the reaction mechanism ---")
    print("The substrate is a secondary alkyl halide and the reagent is a strong, bulky base.")
    print("This combination strongly favors an E2 (bimolecular elimination) reaction.\n")
    
    print("--- Step 2: Analyze the stereochemistry of the substrate ---")
    print(f"The stereochemistry of the substrate is (1S,2R).")
    print("For a 1,2-disubstituted cyclohexane, an (S,R) or (R,S) configuration indicates a 'cis' isomer.")
    print("Therefore, the starting material is cis-1-bromo-2-methylcyclohexane.\n")

    print("--- Step 3: Identify the reactive chair conformation ---")
    print("The E2 reaction requires the leaving group (Br) and a beta-hydrogen to be anti-periplanar (i.e., both axial).")
    print("For a cis-1,2 isomer, this is only possible in the conformation where Br is axial and the methyl group is equatorial.\n")
    
    print("--- Step 4: Identify possible elimination pathways ---")
    print("In the reactive conformation (Br-axial, Me-equatorial), we check the adjacent carbons (C2 and C6) for axial hydrogens:")
    print("  - At C2: The hydrogen is axial. Elimination here forms a double bond between C1 and C2.")
    print("  - At C6: The hydrogen is axial. Elimination here forms a double bond between C1 and C6.\n")

    print("--- Step 5: Determine the possible products ---")
    zaitsev_product = "1-methylcyclohexene"
    hofmann_product = "3-methylcyclohexene"
    print(f"  - Pathway 1 (via C2-H) gives '{zaitsev_product}', the more substituted Zaitsev product.")
    print(f"  - Pathway 2 (via C6-H) gives '{hofmann_product}', the less substituted Hofmann product.\n")

    print("--- Step 6: Apply the selectivity rule for the bulky base ---")
    print(f"The reagent, {reagent}, is a bulky base. Bulky bases favor the formation of the less substituted (Hofmann) product")
    print("by abstracting the more sterically accessible proton on C6.\n")

    # Final Conclusion
    final_product = hofmann_product
    print("### Conclusion ###")
    print(f"The major product of the reaction is {final_product}.")


if __name__ == "__main__":
    identify_elimination_product()