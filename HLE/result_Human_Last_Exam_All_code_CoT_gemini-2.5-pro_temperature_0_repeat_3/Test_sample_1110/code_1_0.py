import textwrap

def explain_reaction_issue():
    """
    Explains the chemical problem and identifies the best solution.
    """
    print("### Analysis of the Synthesis of (2-bromo-4-chlorophenyl)boronic acid ###")
    print("-" * 70)

    # Define the aryl group for clarity in the equations
    aryl_group = "(2-Br-4-Cl-C6H3)"
    print(f"Let's denote the aryl group as Ar = {aryl_group}\n")

    # --- Main Reaction ---
    print("1. The Desired Reaction Pathway (Forms Product 1)")
    print("This pathway should yield a single product with one Boron NMR signal.")
    print("\nEquation 1: Formation of the aryllithium intermediate via halogen-metal exchange.")
    print(f"   1 Ar-I + 1 n-BuLi  ->  1 Ar-Li + 1 n-BuI")
    
    print("\nEquation 2: Reaction with trimethyl borate and subsequent workup.")
    print(f"   1 Ar-Li + 1 B(OMe)3  ->  1 [Ar-B(OMe)3]-Li+")
    print(f"   1 [Ar-B(OMe)3]-Li+ + H+/H2O  ->  1 Ar-B(OH)2  + 3 MeOH + Li+")
    print(f"   Desired Product (Product 1): {aryl_group}-B(OH)2")
    print("-" * 70)

    # --- Side Reaction ---
    print("2. The Problem: A Second Boron Signal (Formation of Product 2)")
    print("The second B signal implies an undesired boron-containing byproduct is also formed.")
    
    explanation = textwrap.fill(
        "This byproduct is most likely the diarylborinic acid, which forms when the highly reactive aryllithium (Ar-Li) attacks the desired boronic ester product instead of the trimethyl borate reagent.",
        width=70
    )
    print(explanation)

    print("\nEquation 3: The Side Reaction Pathway")
    print(f"   1 Ar-Li + 1 Ar-B(OH)2  ->  [ (Ar)2-B(OH)2 ]-Li+")
    print(f"   [ (Ar)2-B(OH)2 ]-Li+ + H+/H2O -> 1 (Ar)2-B(OH) + LiOH + H2")
    print(f"   Byproduct (Product 2): ({aryl_group})2-B(OH)")
    print("-" * 70)

    # --- Analysis and Solution ---
    print("3. Diagnosis and Solution")
    diagnosis = textwrap.fill(
        "This side reaction occurs because of the 1.05 equivalents of n-BuLi used. This slight excess means that after all the starting material (Ar-I) is consumed, there is still some reactive Ar-Li left. This leftover Ar-Li then attacks the product, creating the byproduct.",
        width=70
    )
    print(diagnosis)

    print("\nEvaluating the choices:")
    print(" - A (Decrease temperature): The reaction is already at -78°C; further cooling is unlikely to solve a stoichiometry problem.")
    print(" - B (Use triethylborate): Might slightly reduce the side reaction due to steric hindrance, but doesn't fix the root cause.")
    print(" - D (Use less trimethyl borate): This would make the problem WORSE by making the desired reagent more scarce.")
    print(" - E (Change solvent): THF is the standard, correct solvent for this reaction.")
    print("\n - C (Use more precise amount of n-BuLi): This is the correct solution. By using exactly 1.00 equivalent (or even slightly less), you ensure that the aryllithium (Ar-Li) is the limiting reagent and gets completely consumed by the large excess of trimethyl borate, preventing the formation of the diaryl byproduct.")

explain_reaction_issue()