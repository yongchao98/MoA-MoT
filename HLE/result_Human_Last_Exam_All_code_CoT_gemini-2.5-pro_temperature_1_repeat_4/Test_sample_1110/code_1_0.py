def analyze_synthesis(tmb_eq):
    """
    Models the stoichiometry of the boronic acid synthesis to identify the
    source of the second boron NMR signal.
    
    Args:
        tmb_eq (float): The number of equivalents of trimethyl borate used.
    """
    # For calculation purposes, let's assume we start with 1.0 mole of the starting material.
    starting_material_moles = 1.0
    nBuLi_moles = 1.05  # 1.05 equivalents are used

    # Step 1: Lithiation. The most reactive C-I bond reacts with n-BuLi.
    # Equation: Ar-I + n-BuLi -> Ar-Li + n-BuI
    # Stoichiometry is 1:1. Since n-BuLi is in slight excess, all Ar-I reacts.
    aryllithium_moles_formed = starting_material_moles
    
    # Step 2: Borylation. The aryllithium reacts with trimethyl borate (TMB).
    # Equation: Ar-Li + B(OMe)3 -> [Ar-B(OMe)3]-Li+
    # Stoichiometry is 1:1.
    tmb_initial_moles = tmb_eq * starting_material_moles
    tmb_reacted_moles = aryllithium_moles_formed
    tmb_excess_moles = tmb_initial_moles - tmb_reacted_moles
    
    # Step 3: Aqueous Workup (Hydrolysis).
    # The boronate ester intermediate forms the desired product.
    # Equation: [Ar-B(OMe)3]-Li+ --(H2O)--> Ar-B(OH)2
    product_moles = tmb_reacted_moles
    
    # The excess, unreacted TMB hydrolyzes to boric acid.
    # Equation: B(OMe)3 + 3 H2O -> B(OH)3 + 3 MeOH
    byproduct_moles = tmb_excess_moles
    
    print(f"\n--- Stoichiometric Analysis for {tmb_eq:.1f} eq. of Trimethyl Borate ---")
    print(f"Based on the reaction sequence, the final boron-containing species are:")
    print(f"  - Moles of Desired Product (Ar-B(OH)2): {product_moles:.2f}")
    print(f"  - Moles of Byproduct (Boric Acid, B(OH)3): {byproduct_moles:.2f}")
    print("\nDiagnosis:")
    print("The final crude mixture contains two distinct boron compounds: the desired boronic acid and the byproduct boric acid.")
    print("These two species will give two different signals in the 11B NMR spectrum, which explains the observation.")
    if tmb_eq > 2:
        print("Because a large excess of trimethyl borate was used, the boric acid byproduct is present in a very high amount.")

print("This script analyzes the chemical synthesis to find the cause of the observed problem (two Boron signals in NMR).")

# --- Case 1: Analysis of the situation described in the problem ---
analyze_synthesis(tmb_eq=5.0)

print("\n\n--- Proposed Solution and Re-evaluation ---")
print("The problem is the large quantity of boric acid byproduct formed from the excessive use of trimethyl borate.")
print("The solution is to use a much smaller, yet sufficient, excess of the reagent.")
print("Let's re-evaluate the outcome using a more standard 1.2 equivalents of trimethyl borate.")

# --- Case 2: Analysis using the proposed solution ---
analyze_synthesis(tmb_eq=1.2)

print("\n\n--- Final Conclusion ---")
print("By reducing the trimethyl borate from 5.0 to 1.2 equivalents, the molar ratio of byproduct-to-product decreases from 4.0:1.0 to 0.2:1.0.")
print("This dramatically improves the purity of the crude product and will resolve the issue of having a large, secondary boron signal in the NMR.")
print("Therefore, the most effective action is to use less trimethyl borate.")
