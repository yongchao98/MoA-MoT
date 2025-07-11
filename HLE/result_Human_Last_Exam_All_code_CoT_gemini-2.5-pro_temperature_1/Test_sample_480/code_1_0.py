def solve_maccartney_inference():
    """
    Calculates the final projected natural logic operator for the given
    inference using MacCartney's compositional framework.
    """
    # Define the 7 semantic relations
    relations = {
        'e': {'symbol': '≡', 'name': 'equivalence'},
        'f': {'symbol': '⊏', 'name': 'forward entailment'},
        'r': {'symbol': '⊐', 'name': 'reverse entailment'},
        'a': {'symbol': '^', 'name': 'alternation'},
        'i': {'symbol': '|', 'name': 'independence'},
        'c': {'symbol': '#', 'name': 'cover'},
        'n': {'symbol': '~', 'name': 'negation'},
    }
    # Using a corrected join table based on subsequent research.
    # Keys are row, then column, based on the symbols: ≡, ⊏, ⊐, ^, |, #, ~
    # Short symbols: e, f, r, a, i, c, n
    join_table = {
        # R1 = ≡
        'e': {'e': 'e', 'f': 'f', 'r': 'r', 'a': 'a', 'i': 'i', 'c': 'c', 'n': 'n'},
        # R1 = ⊏
        'f': {'e': 'f', 'f': 'f', 'r': 'i', 'a': 'f', 'i': 'i', 'c': 'c', 'n': 'n'},
        # R1 = ⊐
        'r': {'e': 'r', 'f': 'i', 'r': 'r', 'a': 'a', 'i': 'i', 'c': 'r', 'n': 'n'},
        # R1 = ^
        'a': {'e': 'a', 'f': 'r', 'r': 'a', 'a': 'a', 'i': 'r', 'c': 'a', 'n': 'c'},
        # R1 = |
        'i': {'e': 'i', 'f': 'i', 'r': 'i', 'a': 'i', 'i': 'i', 'c': 'i', 'n': 'c'},
        # R1 = #
        'c': {'e': 'c', 'f': 'c', 'r': 'r', 'a': 'a', 'i': 'i', 'c': 'c', 'n': 'n'},
        # R1 = ~
        'n': {'e': 'n', 'f': 'r', 'r': 'f', 'a': 'c', 'i': 'c', 'c': 'i', 'n': 'e'},
    }

    # Edits from Premise to Hypothesis and their atomic relations
    edits = [
        {'desc': 'Insertion of "not"', 'relation_key': 'n'},          # Negation (~)
        {'desc': 'Deletion of "pop"', 'relation_key': 'f'},           # Forward Entailment (⊏)
        {'desc': 'Substitution of "Taylor Swift" with "Michael Jackson"', 'relation_key': 'a'} # Alternation (^)
    ]

    # Start with the identity relation
    current_relation_key = 'e'
    print(f"Initial Relation: {relations[current_relation_key]['name']} ({relations[current_relation_key]['symbol']})")
    print("-" * 20)

    for i, edit in enumerate(edits):
        context_relation_key = current_relation_key
        edit_relation_key = edit['relation_key']

        # Perform the composition using the join table
        current_relation_key = join_table[context_relation_key][edit_relation_key]

        print(f"Step {i+1}: {edit['desc']}")
        print(f"  Context Relation : {relations[context_relation_key]['name']} ({relations[context_relation_key]['symbol']})")
        print(f"  Edit Relation    : {relations[edit_relation_key]['name']} ({relations[edit_relation_key]['symbol']})")
        print(f"  Projected Relation = join({relations[context_relation_key]['symbol']}, {relations[edit_relation_key]['symbol']}) = {relations[current_relation_key]['symbol']}")
        print(f"  Resulting Relation: {relations[current_relation_key]['name']} ({relations[current_relation_key]['symbol']})")
        print("-" * 20)

    final_relation_name = relations[current_relation_key]['name']
    print(f"The final projected natural logic operator is: {final_relation_name}")

solve_maccartney_inference()