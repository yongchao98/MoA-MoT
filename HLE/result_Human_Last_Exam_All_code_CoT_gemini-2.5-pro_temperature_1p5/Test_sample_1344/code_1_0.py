import textwrap

def suggest_catalyst_system():
    """
    Outlines a promising single-site catalyst system capable of both
    olefin polymerization and polyolefin hydrogenolysis.
    """
    catalyst_system = {
        "Concept": "Switchable Dual-Function Single-Site Catalyst",
        "Metal Center": {
            "Element": "Zirconium (Zr)",
            "Reason": "Excellent activity for both C-C bond formation (polymerization) and C-C bond cleavage (hydrogenolysis)."
        },
        "Ligand_System": {
            "Type": "Surface-Grafted Organometallic Complex (Ligand-Free Approach)",
            "Precursor": "Tetrakis(neopentyl)zirconium, Zr(CH2C(CH3)3)4",
            "Description": "The silica surface itself acts as a bulky, immobilizing ligand. The neopentyl groups are ancillary ligands that are replaced or activated during the catalytic cycles."
        },
        "Support": {
            "Material": "High Surface Area Silica (SiO2)",
            "Pre-treatment": "Heated to 500 degrees Celsius under vacuum.",
            "Purpose": "Controls the density of surface silanol (Si-OH) groups to create isolated, well-defined single-sites via grafting. Provides stability and prevents catalyst dimerization."
        },
        "Active Site Formation": {
            "Reaction": "Zr(CH2C(CH3)3)4 + ≡Si-OH -> (≡SiO)Zr(CH2C(CH3)3)3 + CH4C(CH3)3",
            "Result": "A covalently tethered, uniform single-site zirconium catalyst on the silica surface."
        },
        "Functional Mode 1: Polymerization": {
            "Conditions": "Lower temperatures (e.g., 70-90 °C) with a continuous feed of olefin monomer (e.g., ethylene).",
            "Mechanism": "The grafted zirconium-alkyl site acts as an initiator for classical migratory insertion of olefins, leading to the growth of long polymer chains."
        },
        "Functional Mode 2: Hydrogenolysis (Plastic Breakdown)": {
            "Conditions": "Higher temperatures (e.g., 150-300 °C) under hydrogen (H2) pressure.",
            "Mechanism": "The catalyst forms a zirconium hydride species. This hydride actively breaks down polyolefin chains into short, saturated alkanes (e.g., propane, butane, etc.) through a sequence of C-H activation, β-alkyl elimination, and hydrogenation steps."
        }
    }

    print("--- Proposed Optimal Bifunctional Catalyst System ---")
    for key, value in catalyst_system.items():
        print(f"\n{key.replace('_', ' ')}:")
        if isinstance(value, dict):
            for sub_key, sub_value in value.items():
                wrapped_text = textwrap.fill(f"{sub_key}: {sub_value}", width=80, initial_indent="  ", subsequent_indent="    ")
                print(wrapped_text)
        else:
            print(f"  {value}")

suggest_catalyst_system()
<<<A single-site zirconium catalyst, prepared by grafting an organozirconium complex onto a silica support pre-treated at 500 °C, capable of switching function between polymerization (at <100°C with olefins) and hydrogenolysis (at >150°C with H2).>>>