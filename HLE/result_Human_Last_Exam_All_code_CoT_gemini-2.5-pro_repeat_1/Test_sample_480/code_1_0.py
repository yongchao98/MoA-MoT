def solve_maccartney_inference():
    """
    Solves the natural language inference problem using MacCartney's framework,
    printing each step of the calculation.
    """
    # MacCartney's 7 semantic relations and their symbols
    RELATIONS = {
        '=': 'Equivalence',
        '<': 'Forward Entailment',
        '>': 'Reverse Entailment',
        '^': 'Cover',
        'v': 'Exhaustivity',
        '|': 'Alternation',
        '#': 'Independence'
    }

    # MacCartney's composition (join) table as a dictionary
    # Rows: R1 (accumulated), Columns: R2 (new)
    JOIN_TABLE = {
        #      =    <    >    ^    v    |    #
        '=': ['=', '<', '>', '^', 'v', '|', '#'],
        '<': ['<', '<', '#', '^', '<', '|', '#'],
        '>': ['>', '#', '>', '>', 'v', '>', '#'],
        '^': ['^', '<', '>', '^', '<', '|', '<'],
        'v': ['v', 'v', '>', '>', 'v', '>', '>'],
        '|': ['|', 'v', '<', '<', 'v', '|', '<'],
        '#': ['#', '#', '#', '<', '>', '<', '#'],
    }

    # Relation flipping for downward-monotone contexts
    FLIP_TABLE = {
        '=': '=', '<': '>', '>': '<', '^': 'v', 'v': '^', '|': '|', '#': '#'
    }

    def join(r1, r2):
        """Composes two relations using the join table."""
        cols = ['=', '<', '>', '^', 'v', '|', '#']
        col_idx = cols.index(r2)
        return JOIN_TABLE[r1][col_idx]

    def flip(r):
        """Flips a relation for a downward-monotone context."""
        return FLIP_TABLE[r]

    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    print("Analyzing the inference:")
    print(f"Premise (P): \"{premise}\"")
    print(f"Hypothesis (H): \"{hypothesis}\"")
    print("-" * 20)
    print("Applying edits from left to right to transform P to H:\n")

    # Step 0: Initialization
    accumulated_relation = '='
    print(f"Step 0: Start with the premise itself. The initial relation is self-entailment, which is Equivalence ('=').")
    print(f"Accumulated Relation: {accumulated_relation}\n")

    # Step 1: First edit ("is singing" -> "is not singing")
    print("Step 1: Process the first edit: \"is singing\" -> \"is not singing\".")
    atomic_rel1 = '|'
    print(f"  - The atomic relation is Alternation ('{atomic_rel1}'), as 'singing' and 'not singing' are mutually exclusive.")
    # The context is neutral (upward-monotone), so no flip is needed.
    projected_rel1 = atomic_rel1
    print(f"  - The context is upward-monotone, so the projected relation is unchanged: '{projected_rel1}'.")
    
    prev_relation = accumulated_relation
    accumulated_relation = join(prev_relation, projected_rel1)
    print(f"  - Compose the relations: join({prev_relation}, {projected_rel1}) = {accumulated_relation}")
    print(f"Accumulated Relation is now '{accumulated_relation}' ({RELATIONS[accumulated_relation]}).\n")

    # Step 2: Second edit ("a pop song by Taylor Swift" -> "a song by Michael Jackson")
    print("Step 2: Process the second edit: \"... a pop song by Taylor Swift\" -> \"... a song by Michael Jackson\".")
    atomic_rel2 = '|'
    print(f"  - The atomic relation is Alternation ('{atomic_rel2}'), as the sets of songs by these two artists are disjoint.")
    
    # The context is "Mark is not singing [...]", which is downward-monotone due to "not".
    projected_rel2 = flip(atomic_rel2)
    print(f"  - This edit occurs in a downward-monotone context ('not...'). We must flip the atomic relation: flip({atomic_rel2}) = {projected_rel2} ({RELATIONS[projected_rel2]}).")

    prev_relation = accumulated_relation
    accumulated_relation = join(prev_relation, projected_rel2)
    print(f"  - Compose the accumulated relation with this new projected relation.")
    print(f"  - Final Equation: join({prev_relation}, {projected_rel2}) = {accumulated_relation}")
    print("-" * 20)
    
    final_relation_name = RELATIONS[accumulated_relation]
    print(f"The final projected natural logic operator is '{accumulated_relation}'.")
    print(f"The name of this operator is: {final_relation_name}")

if __name__ == '__main__':
    solve_maccartney_inference()