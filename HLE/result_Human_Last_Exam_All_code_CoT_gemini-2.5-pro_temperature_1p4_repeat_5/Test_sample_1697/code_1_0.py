import textwrap

def identify_reaction_product():
    """
    This script determines the product of a two-step organic reaction
    using chemical principles.
    """
    print("--- Task: Identify the chemical compound ---")
    print("Reaction: N,N-diethyl-3-dimethylaminobenzamide reacts with 1) sec-BuLi/TMEDA in THF, then 2) methyl iodide.\n")

    print("--- Step-by-Step Analysis ---")

    # Step 1: Directed ortho-Metalation (DoM)
    print("\n[Step 1: Lithiation with sec-BuLi and TMEDA]")
    print("Starting Material: N,N-diethyl-3-dimethylaminobenzamide")
    print("This molecule has two key functional groups on the benzene ring:")
    print("  - A N,N-diethylamide group [-C(=O)N(Et)2] at position 1.")
    print("  - A dimethylamino group [-N(Me)2] at position 3.")
    
    print("\nReagents: sec-BuLi is a strong base, and TMEDA is a chelating agent.")
    
    explanation_step1 = """
    These reagents perform a 'Directed ortho-Metalation' (DoM), where a proton
    on the ring is removed. The position of removal is directed by the functional
    groups, which are called Directing Metalation Groups (DMGs).

    - The amide group at position 1 directs the base to remove a proton from either position 2 or 6.
    - The amino group at position 3 directs the base to remove a proton from either position 2 or 4.

    Both groups direct the reaction to position 2. This 'cooperative effect' makes the
    proton at position 2 the most acidic. Therefore, the sec-BuLi selectively
    removes the proton at position 2, replacing it with a lithium atom.
    """
    print(textwrap.dedent(explanation_step1))
    print("Intermediate formed: 2-lithio-N,N-diethyl-3-dimethylaminobenzamide.\n")

    # Step 2: Electrophilic Quench
    print("[Step 2: Reaction with Methyl Iodide]")
    print("The aryllithium intermediate formed in Step 1 is a strong nucleophile (a carbanion at position 2).")
    print("Methyl iodide (CH3I) is an electrophile.")

    explanation_step2 = """
    The nucleophilic carbon at position 2 attacks the electrophilic methyl group of
    methyl iodide. This forms a new carbon-carbon bond, adding a methyl group
    to position 2 of the benzene ring.
    """
    print(textwrap.dedent(explanation_step2))

    # Final Product
    print("--- Conclusion ---")
    final_product_name = "N,N-diethyl-2-methyl-3-dimethylaminobenzamide"
    print(f"The final product is the starting material with a methyl group added at position 2.")
    print(f"\nFinal Product Name: {final_product_name}\n")
    
    print("Numbers in the final chemical name:")
    print("Position number for the new methyl group: 2")
    print("Position number for the existing dimethylamino group: 3")
    print("The designation 'N,N' indicates that the two ethyl groups are attached to the Nitrogen atom of the amide.")

identify_reaction_product()
<<<N,N-diethyl-2-methyl-3-dimethylaminobenzamide>>>