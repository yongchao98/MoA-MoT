def solve_entailment():
    """
    Calculates and explains the projected natural logic operator for the given inference
    based on MacCartney's framework.
    """
    
    # Define the 7 logic relations and their symbols
    relations = {
        '≡': "Equivalence",
        '<': "Forward Entailment",
        '>': "Reverse Entailment",
        '|': "Alternation",
        '⌣': "Cover",
        '^': "Negation",
        '#': "Independence"
    }

    # B(R1, R2) composition table, including only the entries needed for this problem.
    # The table defines the result of composing relation R1 with R2 (R1 ◦ R2).
    composition_table = {
        ('≡', '^'): '^',  # Equivalence ◦ Negation = Negation
        ('^', '>'): '|',  # Negation ◦ Reverse Entailment = Alternation
        ('|', '⌣'): '<'   # Alternation ◦ Cover = Forward Entailment
    }

    # op(f, c) table for effect of context 'c' on a lexical relation 'f'.
    # '¬' represents a downward-monotone (negative) context.
    op_table = {
        ('<', '¬'): '>',  # Deletion in a negative context inverts forward entailment to reverse.
        ('|', '¬'): '⌣'   # Substitution of disjoint items in a negative context turns alternation to cover.
    }

    # The sequence of edits transforming premise to hypothesis
    steps = [
        {
            "description": "Insertion of 'not'",
            "base_relation": '^',
            "context": None,
        },
        {
            "description": "Deletion of 'pop'",
            "base_relation": '<',
            "context": '¬',
        },
        {
            "description": "Substitution of 'Taylor Swift' with 'Michael Jackson'",
            "base_relation": '|',
            "context": '¬',
        }
    ]

    print("Deriving the entailment relation step-by-step:\n")
    print("Premise: 'Mark is singing a pop song by Taylor Swift'")
    print("Hypothesis: 'Mark is not singing a song by Michael Jackson'\n")
    
    # Start with the Equivalence relation (premise vs. itself)
    current_relation = '≡'
    print(f"Initial state: Relation is {relations[current_relation]} ({current_relation})\n" + "-"*40)

    for i, step in enumerate(steps, 1):
        # Determine the relation for the edit, considering context
        if step["context"]:
            edit_relation = op_table[(step["base_relation"], step["context"])]
        else:
            edit_relation = step["base_relation"]
        
        # Compose the current relation with the new edit's relation
        new_relation = composition_table[(current_relation, edit_relation)]

        print(f"Step {i}: {step['description']}")
        
        # Show the symbolic composition equation
        print(f"  Equation: {relations[current_relation]} ({current_relation}) ◦ {relations[edit_relation]} ({edit_relation}) = {relations[new_relation]} ({new_relation})")
        
        current_relation = new_relation
        print(f"  New projected relation: {relations[current_relation]}\n" + "-"*40)

    final_relation_name = relations[current_relation]
    print(f"\nThe final projected natural logic operator is: {final_relation_name}")


solve_entailment()