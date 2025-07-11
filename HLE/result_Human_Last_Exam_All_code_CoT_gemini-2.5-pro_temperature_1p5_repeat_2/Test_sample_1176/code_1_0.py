def solve_geology_quiz():
    """
    Analyzes geological statements and prints the result as a string.

    Each statement is classified as Consensus (C) or Debated (D).
    The function prints the classification for each statement and then
    the final concatenated string of answers.
    """
    # Statements and their classification (C=Consensus, D=Debated)
    classifications = {
        1: ("C", "The Morrison formation represents a foredeep deposit."),
        2: ("D", "Metamorphic core complexes formed in response to a slab window."),
        3: ("C", "An orogenic plateau ('Nevadaplano') existed by the end of the Cretaceous."),
        4: ("D", "An 'Arizonaplano' was constructed via magmatic thickening."),
        5: ("C", "Laramide-style structures are generally inboard of Sevier-style structures."),
        6: ("C", "The Sierra Nevada and Idaho Batholiths were formed by the Cordilleran arc."),
        7: ("C", "A series of large-magnitude ignimbrite eruptions propagated to the southwest."),
        8: ("D", "The Colorado Plateau reached modern elevations by the late Cretaceous."),
        9: ("D", "The Cascades arc reached modern elevation by the late Oligocene."),
        10: ("D", "Laramide-style deformation was caused by subduction of the Shatsky conjugate.")
    }

    final_answer_string = ""

    print("Analysis of each statement:")
    # Loop through the statements in order and print the breakdown
    for i in sorted(classifications.keys()):
        code, _ = classifications[i]
        # This step prints each number and its corresponding letter,
        # similar to the elements of an equation.
        print(f"Statement {i}: {code}")
        final_answer_string += code
    
    print("\n" + "="*30)
    print("Final Result String:")
    print(final_answer_string)
    print("="*30)

    # Final answer in the required format
    print(f"<<<{final_answer_string}>>>")

solve_geology_quiz()