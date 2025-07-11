import textwrap

def present_optimal_catalyst_system():
    """
    Presents a state-of-the-art catalyst system for both olefin polymerization
    and polyolefin hydrogenolysis.
    """
    
    metal = {
        "name": "Zirconium (Zr)",
        "group": "Group IV",
        "rationale": (
            "Zirconium is a classic Group IV metal renowned for its high activity "
            "in olefin polymerization (as seen in Zirconocene catalysts). "
            "Recent research has shown that when properly supported, it can also "
            "catalyze the cleavage of C-C bonds required for hydrogenolysis."
        )
    }

    ligand_system = {
        "name": "Surface Siloxy (Si-O-) and Hydride (H-) Ligands",
        "type": "Single-Site, Post-Metallocene",
        "rationale": (
            "Instead of a traditional cyclopentadienyl ligand, the active site is "
            "formed by grafting a simple zirconium precursor (like Zr(CH2tBu)4) onto "
            "a silica surface. The actual active site becomes a surface-bound "
            "zirconium hydride, (≡SiO)x-Zr-H. This site is highly reactive and "
            "undercoordinated, allowing it to perform both monomer insertion for "
            "polymerization and interaction with the polymer backbone for cleavage."
        )
    }

    support = {
        "name": "Silica (SiO2)",
        "type": "Inorganic Oxide Support",
        "rationale": (
            "The silica support is not inert. It is essential for creating and "
            "stabilizing the highly reactive single-site zirconium hydride species. "
            "It prevents the catalyst from deactivating via dimerization and provides "
            "the necessary structural and electronic environment for the dual "
            "catalytic function."
        )
    }
    
    print("="*60)
    print("Proposed Dual-Function Catalyst System")
    print("="*60)
    
    print(f"\n[ Metal ]")
    print(f"Choice: {metal['name']} ({metal['group']})")
    print("\nRationale:")
    print(textwrap.fill(metal['rationale'], width=60))
    
    print("\n" + "-"*60 + "\n")

    print(f"[ Ligand System ]")
    print(f"Choice: {ligand_system['name']}")
    print(f"Type: {ligand_system['type']}")
    print("\nRationale:")
    print(textwrap.fill(ligand_system['rationale'], width=60))

    print("\n" + "-"*60 + "\n")
    
    print(f"[ Support ]")
    print(f"Choice: {support['name']}")
    print(f"Type: {support['type']}")
    print("\nRationale:")
    print(textwrap.fill(support['rationale'], width=60))

    print("\n" + "="*60)
    print("Summary:")
    print("The optimal combination is a single-site zirconium hydride")
    print("catalyst supported on silica (e.g., ZrH/SiO2).")
    print("Reaction control (polymerization vs. hydrogenolysis) can be")
    print("tuned by varying temperature, and pressures of H2 and olefin.")
    print("="*60)


if __name__ == '__main__':
    present_optimal_catalyst_system()
