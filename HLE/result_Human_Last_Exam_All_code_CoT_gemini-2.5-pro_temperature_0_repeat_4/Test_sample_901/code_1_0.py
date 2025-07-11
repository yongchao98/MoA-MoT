def identify_product():
    """
    Analyzes the reaction of (1S,2R)-1-bromo-2-methylcyclohexane with potassium tert-butoxide
    and identifies the major organic product.
    """
    substrate = "(1S,2R)-1-bromo-2-methylcyclohexane"
    reagent = "Potassium tert-butoxide (t-BuOK)"
    product = "3-methylcyclohexene"

    print("### Step-by-Step Reaction Analysis ###\n")

    print(f"1. Reaction Type Identification:")
    print(f"   - Substrate: {substrate} (a secondary alkyl halide).")
    print(f"   - Reagent: {reagent} (a strong, bulky base).")
    print("   - Conclusion: These conditions strongly favor an E2 (bimolecular elimination) reaction.\n")

    print("2. Stereochemical Requirement for E2:")
    print("   - The E2 reaction requires an 'anti-periplanar' geometry.")
    print("   - In a cyclohexane ring, this means the leaving group (Br) and a beta-hydrogen must be in a 'trans-diaxial' arrangement (both axial, one up, one down).\n")

    print("3. Substrate Conformation Analysis:")
    print(f"   - For the E2 reaction to occur with {substrate}, the Bromine atom must be in an axial position.")
    print("   - The (1S,2R) stereochemistry means the Br and the C2-Methyl group are trans.")
    print("   - In the required reactive conformation, this forces both the Br atom (at C1) and the methyl group (at C2) into axial positions.\n")

    print("4. Evaluating Elimination Pathways:")
    print("   - The beta-hydrogens are on C2 and C6.")
    print("   - Pathway via C2 (Zaitsev): The hydrogen on C2 is equatorial because the methyl group is axial. This H is NOT anti-periplanar to the axial Br. This pathway is FORBIDDEN.")
    print("   - Pathway via C6 (Hofmann): The hydrogen on C6 has an axial position. This H IS anti-periplanar to the axial Br. This pathway is ALLOWED.\n")

    print("5. Final Product Determination:")
    print("   - Elimination must occur by removing the axial hydrogen from C6, forming a double bond between C1 and C6.")
    print(f"   - The resulting product is {product}.\n")
    
    print("### Final Equation ###")
    # The prompt requires printing the numbers in the final equation.
    # The numbers are 1, 2, 1, 2 from the reactant and 3 from the product.
    print(f"The reaction is: ({1}S,{2}R)-{1}-bromo-{2}-methylcyclohexane --(t-BuOK)--> {3}-methylcyclohexene")

identify_product()