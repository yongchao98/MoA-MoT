def analyze_complex_stability():
    """
    Analyzes the stability of four Iridium(III) complexes based on their ligand structure.

    The stability, and thus the operational lifetime in a Light-Emitting Electrochemical Cell (LEC),
    of these cyclometalated Iridium complexes is highly dependent on the substituents on the
    phenylpyridine (C^N) ligands.

    A key principle is that a fluorine atom at the ortho-position (position 2) of the phenyl ring
    (relative to the Ir-C bond) introduces a significant degradation pathway, leading to a shorter
    lifetime. Fluorine atoms at other positions (meta or para) generally enhance stability.

    This function identifies the complexes with this destabilizing feature.
    """

    # Data for the C^N ligand of each complex.
    # 'F_positions' lists the positions of fluorine atoms on the phenyl ring.
    complexes = {
        1: {'name': 'Complex 1', 'ligand_type': '2,4-difluorophenyl', 'F_positions': [2, 4]},
        2: {'name': 'Complex 2', 'ligand_type': 'phenyl', 'F_positions': []},
        3: {'name': 'Complex 3', 'ligand_type': '2-fluorophenyl', 'F_positions': [2]},
        4: {'name': 'Complex 4', 'ligand_type': '3,5-difluorophenyl', 'F_positions': [3, 5]}
    }

    shorter_lifetime_complexes = []
    ortho_position = 2

    print("Analysis of Complex Stability and Expected Lifetime:")
    print("=" * 50)
    print(f"Rule: Complexes with a fluorine atom at the ortho-position ({ortho_position}) are expected to have shorter lifetimes due to known degradation pathways.")
    print("=" * 50)

    for cid, data in complexes.items():
        # Check if the ortho-position is in the list of fluorine positions
        if ortho_position in data['F_positions']:
            status = "Shorter lifetime expected."
            shorter_lifetime_complexes.append(cid)
        else:
            status = "Longer lifetime expected."
        
        print(f"Analyzing {data['name']} (ligand: {data['ligand_type']}-pyridine):")
        print(f"  - Fluorine positions: {data['F_positions'] if data['F_positions'] else 'None'}")
        print(f"  - Conclusion: {status}\n")

    shorter_lifetime_complexes.sort()
    
    print("Final Result:")
    print("The complexes expected to show shorter lifetimes are those with an ortho-fluorine substituent.")
    print("These are complexes:")
    for num in shorter_lifetime_complexes:
        print(num)

# Execute the analysis
analyze_complex_stability()