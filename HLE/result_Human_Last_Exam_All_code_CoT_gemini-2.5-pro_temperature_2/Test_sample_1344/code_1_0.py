import textwrap

def find_optimal_catalyst_system():
    """
    This function analyzes the requirements for a dual-purpose catalyst
    and presents the optimal system based on established chemical principles.
    The catalyst must perform both olefin polymerization and polyolefin hydrogenolysis.
    """

    # --- Component Definition ---
    # Based on chemical knowledge, we define the optimal components.
    metal = "Zirconium (Zr)"
    support = "Silica (SiO2), specifically highly dehydroxylated"
    precursor_ligands = "Alkyl groups, e.g., Neopentyl (-CH2C(CH3)3)"
    active_site_ligand = "Surface Siloxide groups (≡Si-O-)"
    active_species = "Surface-supported Zirconium Hydride"
    
    # --- Output the Final Answer ---
    print("### Optimal Single-Site Catalyst System ###")
    print("-" * 40)
    
    print(f"1. GROUP IV METAL:")
    print(f"   - {metal}")
    rationale_metal = (
        "Zirconium offers an excellent balance of high activity for "
        "olefin polymerization and the ability to activate strong C-H and C-C bonds, "
        "which is essential for hydrogenolysis."
    )
    print(textwrap.fill(rationale_metal, initial_indent='     ', subsequent_indent='     ', width=70))
    print("-" * 40)

    print("2. LIGAND / SUPPORT SYSTEM:")
    print(f"   - Support: {support}")
    print(f"   - Ligand:   {active_site_ligand}")

    rationale_support = (
        "A rigid silica support is crucial. It acts as a robust 'ligand' by "
        "covalently binding the Zirconium atom. This prevents the formation of inactive "
        "dimers and creates a well-defined, isolated 'single-site' that is sterically "
        "accessible and highly reactive."
    )
    print(textwrap.fill(rationale_support, initial_indent='     ', subsequent_indent='     ', width=70))
    print("-" * 40)
    
    print("3. FINAL CATALYTIC SPECIES EQUATION (Formation):")
    # Using descriptive names instead of a formal, balanced equation for clarity.
    print("\n   [Zr(Alkyl)4]  +  [Surface-(Si-OH)]   --->   [Surface-(Si-O)-Zr(Alkyl)3]  +  Alkane")
    print("      (Precursor)    (Silica Surface)               (Grafted Species)")

    print("\n   [Surface-(Si-O)-Zr(Alkyl)3]  +  H2   --->   [Surface-(Si-O)-Zr-H(x)]  +  Alkane")
    print("         (Grafted Species)        (Gas)               (ACTIVE CATALYST)")
    print("-" * 40)

    print("\nCONCLUSION:")
    summary = (
        "The optimal combination is a single-site Zirconium Hydride supported on silica. "
        "This system leverages the support as a key ligand to create a coordinatively "
        "unsaturated and highly reactive metal center that can efficiently catalyze both "
        "the making (polymerization) and breaking (hydrogenolysis) of polyolefins."
    )
    print(textwrap.fill(summary, width=70))


if __name__ == "__main__":
    find_optimal_catalyst_system()