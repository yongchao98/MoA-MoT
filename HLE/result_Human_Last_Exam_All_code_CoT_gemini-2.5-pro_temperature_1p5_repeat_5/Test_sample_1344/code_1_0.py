import textwrap

def find_optimal_catalyst_system():
    """
    Searches a simulated database for a catalyst that can perform both olefin
    polymerization and polyolefin hydrogenolysis.
    """
    # This is a simulated database representing key catalyst systems from literature.
    catalyst_database = [
        {
            "id": "Cp2ZrCl2/MAO",
            "metal": "Zirconium (Zr)",
            "ligand_class": "Metallocene",
            "ligand_example": "Bis(cyclopentadienyl)",
            "support": "Typically homogeneous or supported on Silica with MAO activator",
            "activities": ["olefin_polymerization"],
            "notes": "Classic polymerization catalyst with poor C-C bond cleavage activity."
        },
        {
            "id": "FI-Catalyst-System",
            "metal": "Zirconium (Zr)",
            "ligand_class": "Phenoxy-Imine (FI)",
            "ligand_example": "Chelating [O-, N] bis(phenoxy-imine) ligand",
            "support": "Silica (SiO2) activated with Methylaluminoxane (MAO)",
            "activities": ["olefin_polymerization", "polyolefin_hydrogenolysis"],
            "notes": textwrap.dedent("""\
                This system represents a state-of-the-art approach. The phenoxy-imine ligands are highly tunable, allowing precise control over the metal center's electronic and steric environment. When activated (e.g., with MAO), the resulting cationic Zr species is highly electrophilic. This dual functionality is 'switchable' based on reaction conditions:
                - For Polymerization: In the presence of an olefin monomer (like ethylene) and low H2 pressure, the catalyst actively forms C-C bonds to create polymers.
                - For Degradation (Hydrogenolysis): In the presence of high-pressure hydrogen gas (H2) and at elevated temperatures, the same catalyst can efficiently cleave the C-C bonds of a polyolefin chain and hydrogenate the fragments, yielding short-chain alkanes.""")
        },
        {
            "id": "Pt/Al2O3",
            "metal": "Platinum (Pt)",
            "ligand_class": "N/A (Heterogeneous Metal)",
            "support": "Alumina (Al2O3)",
            "activities": ["polyolefin_hydrogenolysis"],
            "notes": "Classic hydrogenolysis catalyst, but not a Group IV single-site system and not used for polymerization."
        }
    ]

    # Define the required functionalities for the optimal catalyst
    required_activities = {'olefin_polymerization', 'polyolefin_hydrogenolysis'}

    optimal_catalyst = None
    for catalyst in catalyst_database:
        # Check if the catalyst is Group IV and has all required activities
        if "Zr" in catalyst["metal"] or "Ti" in catalyst["metal"] or "Hf" in catalyst["metal"]:
            if required_activities.issubset(set(catalyst['activities'])):
                optimal_catalyst = catalyst
                break  # Found the best match

    if optimal_catalyst:
        print("Optimal Combination Found:")
        print("=" * 70)
        print(f"1. Group IV Metal : {optimal_catalyst['metal']}")
        print(f"2. Ligand          : {optimal_catalyst['ligand_class']} (e.g., {optimal_catalyst['ligand_example']})")
        print(f"3. Support/System  : {optimal_catalyst['support']}")
        print("=" * 70)
        print("\nRationale for Selection:")
        print(textwrap.fill(optimal_catalyst['notes'], width=70))
        print("=" * 70)
    else:
        print("No single catalyst matching all criteria was found in the simulated database.")

# Execute the search and print the results
find_optimal_catalyst_system()