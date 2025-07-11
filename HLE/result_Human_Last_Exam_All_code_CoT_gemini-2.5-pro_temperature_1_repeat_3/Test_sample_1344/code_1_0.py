def propose_catalyst_system():
    """
    Proposes a promising single-site catalyst system for both olefin
    polymerization and polyolefin hydrogenolysis.

    The "optimal" catalyst is a subject of ongoing research. This function
    presents a strong candidate based on current scientific literature,
    combining components known for high performance and stability. The goal is
    a single-site system that can perform both C-C bond formation (polymerization)
    and C-C bond cleavage (hydrogenolysis).
    """

    # --- Component Selection ---

    # 1. Group IV Metal: Zirconium is selected for its high activity.
    metal = "Zirconium (Zr)"

    # 2. Ligand: A robust, post-metallocene ligand is needed to create a stable,
    #    well-defined single active site. Pyridyl-diamido ligands are excellent
    #    candidates.
    ligand = "Pyridyl-diamido Ligand"

    # 3. Support: A support is crucial for practical applications, helping to
    #    stabilize the catalyst and prevent deactivation. Silica is a
    #    versatile and widely used support material.
    support = "Silica (SiO2) Support"

    # --- Output the Proposed Catalyst "Equation" ---
    # The final "equation" is the combination of the selected components.
    # We print each component individually as requested.

    print("Proposed Catalyst System for Dual-Function Polymerization and Hydrogenolysis:")
    print("-" * 70)
    print(f"Component 1 (Group IV Metal): {metal}")
    print(f"Component 2 (Ligand Class):   {ligand}")
    print(f"Component 3 (Support):        {support}")
    print("-" * 70)
    print(f"Final Proposed Combination: {metal} complex with a {ligand} on a {support}")
    print("-" * 70)


if __name__ == "__main__":
    propose_catalyst_system()