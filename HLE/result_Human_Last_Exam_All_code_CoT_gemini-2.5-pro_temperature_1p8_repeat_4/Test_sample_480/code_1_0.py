def solve_mccartney_inference():
    """
    Calculates the final natural logic operator for a given premise and hypothesis
    by composing relations from sequential edits according to MacCartney's framework.
    """
    # Define the 7 operators and their symbols
    RELATIONS = {
        '≡': "Identity",
        '<': "Reverse Entailment (Subsumption)",
        '>': "Forward Entailment",
        '~': "Negation (Contradiction)",
        '|': "Alternation (Disjointness)",
        '^': "Cover (Join/Exhaustiveness)",
        '#': "Independence",
        'v': "Cover (from Stanford materials, sometimes called Exhaustiveness)"
    }
    # MacCartney's 7x7 composition table T[R1, R2] = R1 ○ R2
    # R1 is the row, R2 is the column
    COMPOSITION_TABLE = {
        # R2=    ≡    <    >    ~    |    ^    #
        '≡':   { '≡':'≡', '<':'<', '>':'>', '~':'~', '|':'|', '^':'^', '#':'#' },
        '<':   { '≡':'<', '<':'<', '>':'#', '~':'|', '|':'|', '^':'|', '#':'#' },
        '>':   { '≡':'>', '<':'#', '>':'>', '~':'^', '|':'#', '^':'^', '#':'#' },
        '~':   { '≡':'~', '<':'>', '>':'<', '~':'≡', '|':'^', '^':'|', '#':'#' },
        '|':   { '≡':'|', '<':'|', '>':'v', '~':'^', '|':'v', '^':'v', '#':'v' },
        '^':   { '≡':'^', '<':'v', '>':'|', '~':'v', '|':'v', '^':'v', '#':'v' },
        '#':   { '≡':'#', '<':'#', '>':'#', '~':'#', '|':'#', '^':'#', '#':'#' }
    }
    
    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    print(f"Premise (P): \"{premise}\"")
    print(f"Hypothesis (H): \"{hypothesis}\"\n")
    print("Calculating the final relation by composing relations from left-to-right edits:\n")

    # Edits and their corresponding relations between intermediate steps
    # r_i = Relation(H_{i-1}, H_i)
    edits = [
        ("is singing -> is not singing", '~'),
        ("a pop song -> a song", '<'),
        ("Taylor Swift -> Michael Jackson", '^'),
    ]

    current_relation = '≡'
    
    print(f"Step 0: Initial relation between P and itself.")
    print(f"         R_0 = {current_relation} ({RELATIONS[current_relation]})")
    
    for i, (desc, edit_relation) in enumerate(edits):
        step_num = i + 1
        prev_relation = current_relation
        
        # Composition: R_new = R_old ○ r_edit
        current_relation = COMPOSITION_TABLE[prev_relation][edit_relation]
        
        print(f"\nStep {step_num}: Edit is '{desc}'.")
        print(f"         Relation of this edit step (r_{step_num}) is '{edit_relation}' ({RELATIONS[edit_relation]})")
        print(f"         Composing with previous relation: R_{step_num} = R_{step_num-1} ○ r_{step_num}")
        final_equation = f"         Equation: R_{step_num} = {prev_relation} ○ {edit_relation} = {current_relation}"
        print(final_equation)
        print(f"         New running relation is '{current_relation}' ({RELATIONS[current_relation]})")

    final_operator_symbol = current_relation
    final_operator_name = RELATIONS[final_operator_symbol]

    print("\n---------------------------------------------------------")
    print(f"The final projected natural logic operator is '{final_operator_name}'.")
    print("---------------------------------------------------------")

solve_mccartney_inference()