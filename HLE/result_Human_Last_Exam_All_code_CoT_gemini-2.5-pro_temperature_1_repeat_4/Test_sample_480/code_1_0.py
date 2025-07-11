def solve_natural_logic_inference():
    """
    Calculates the final projected natural logic operator for a given
    premise and hypothesis based on MacCartney's framework.
    """
    # MacCartney's 7 logic relations
    # ≡ (equiv), sq (forward_entail), v (reverse_entail), ! (negation),
    # | (alternation), # (cover), ^ (independence)
    relations = ['≡', 'sq', 'v', '!', '|', '#', '^']
    
    # MacCartney's Join Table: R_new = join(R_current, R_edit)
    # Rows are R_current, Columns are R_edit
    join_table_data = [
        # ≡    sq   v    !    |    #    ^
        ['≡', 'sq', 'v', '!', '|', '#', '^'],  # ≡
        ['sq', 'sq', 'v', '!', 'v', '#', '^'],  # sq
        ['v', 'v', 'v', '|', '|', '#', '^'],  # v
        ['!', 'v', 'sq', '≡', '^', '#', '|'],  # !
        ['|', 'sq', 'sq', '^', '≡', '#', 'v'],  # |
        ['#', '#', '#', '#', '#', '#', '#'],  # #
        ['^', '^', '^', '^', '^', '#', '≡'],  # ^
    ]

    # Create a dictionary-based table for easy lookup
    join_table = {
        outer: {
            inner: result for inner, result in zip(relations, row)
        } for outer, row in zip(relations, join_table_data)
    }

    # Define the sequence of edits and their lexical relations
    edits = [
        ("Negation: 'is singing' -> 'is not singing'", '!'),
        ("Deletion: 'pop song' -> 'song'", 'sq'),
        ("Substitution: 'Taylor Swift' -> 'Michael Jackson'", '|'),
    ]

    # Start with the Identity relation
    current_relation = '≡'
    
    print("Starting with Premise: 'Mark is singing a pop song by Taylor Swift'")
    print(f"Initial relation to itself is Identity ({current_relation})\n")

    # Sequentially apply the edits
    for i, (description, edit_relation) in enumerate(edits):
        previous_relation = current_relation
        current_relation = join_table[previous_relation][edit_relation]
        print(f"Step {i+1}: {description}")
        print(f"  - Lexical relation of edit: {edit_relation}")
        print(f"  - Projecting over current relation '{previous_relation}': join({previous_relation}, {edit_relation}) = {current_relation}")
        print(f"  - New projected relation is: {current_relation}\n")

    relation_names = {
        '≡': 'Identity',
        'sq': 'Forward Entailment',
        'v': 'Reverse Entailment',
        '!': 'Negation',
        '|': 'Alternation',
        '#': 'Cover',
        '^': 'Independence',
    }
    
    final_operator_name = relation_names[current_relation]
    print(f"The final projected natural logic operator is {current_relation}, which is named '{final_operator_name}'.")

solve_natural_logic_inference()