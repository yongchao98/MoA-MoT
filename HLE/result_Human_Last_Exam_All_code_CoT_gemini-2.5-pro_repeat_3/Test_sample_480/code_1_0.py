def solve_maccartney_inference():
    """
    Calculates the final projected natural logic operator for a given
    premise-hypothesis pair based on MacCartney's framework.
    """
    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    # Define the sequence of edits and their atomic relations
    # < : Forward Entailment
    # | : Alternation (mutually exclusive)
    # ^ : Negation
    edits = [
        {'description': 'Deletion of "pop"', 'relation': '<'},
        {'description': 'Substitution: "Taylor Swift" -> "Michael Jackson"', 'relation': '|'},
        {'description': 'Insertion of "not"', 'relation': '^'}
    ]

    # Define relation symbols and their names
    relation_names = {
        '=': 'Equivalence',
        '<': 'Forward Entailment',
        '>': 'Reverse Entailment',
        '^': 'Negation',
        '|': 'Alternation',
        'v': 'Cover',  # Using 'v' for the union symbol ∪
        '#': 'Independence'
    }

    # MacCartney's join table for composing relations
    # Format: join_table[current_relation][edit_relation]
    join_table = {
        '=': {'=': '=', '<': '<', '>': '>', '^': '^', '|': '|', 'v': 'v', '#': '#'},
        '<': {'=': '<', '<': '<', '>': '#', '^': '>', '|': '|', 'v': '<', '#': '#'},
        '>': {'=': '>', '<': '#', '>': '>', '^': '<', '|': '>', 'v': 'v', '#': '#'},
        '^': {'=': '^', '<': 'v', '>': '|', '^': '=', '|': '>', 'v': '<', '#': '#'},
        '|': {'=': '|', '<': '|', '>': 'v', '^': '<', '|': '#', 'v': '|', '#': '#'},
        'v': {'=': 'v', '<': '<', '>': 'v', '^': '>', '|': 'v', 'v': '#', '#': '#'},
        '#': {'=': '#', '<': '#', '>': '#', '^': '#', '|': '#', 'v': '#', '#': '#'},
    }

    # Start with the Equivalence relation
    current_relation = '='

    print(f"Premise: \"{premise}\"")
    print(f"Hypothesis: \"{hypothesis}\"")
    print("-" * 30)
    print(f"Calculating the final projected relation step-by-step:")
    print(f"1. Initial Relation: {current_relation} ({relation_names[current_relation]})")

    # Sequentially compose the relations for each edit
    for i, edit in enumerate(edits):
        old_relation = current_relation
        edit_relation = edit['relation']
        current_relation = join_table[old_relation][edit_relation]
        
        print(f"\nStep {i + 2}: Applying edit '{edit['description']}'")
        print(f"  - Atomic relation of edit is: {edit_relation} ({relation_names[edit_relation]})")
        # Outputting the 'equation' as requested
        print(f"  - Composition equation: join({old_relation}, {edit_relation}) = {current_relation}")
        print(f"  - New projected relation is: {current_relation} ({relation_names[current_relation]})")

    print("-" * 30)
    final_answer_name = relation_names[current_relation]
    print(f"The name of the final projected natural logic operator is: {final_answer_name}")
    print(f"<<<{final_answer_name}>>>")

if __name__ == "__main__":
    solve_maccartney_inference()