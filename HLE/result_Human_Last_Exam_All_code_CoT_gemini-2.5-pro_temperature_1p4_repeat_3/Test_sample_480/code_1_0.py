def solve_maccartney_inference():
    """
    Determines the final projected natural logic operator for a given inference
    using MacCartney's compositional semantics.
    """
    # Define the 7 semantic relations and their names
    RELATIONS = {
        '=': "Equivalence",
        'sq': "Forward Entailment",
        '^': "Reverse Entailment",
        '|': "Alternation",
        'v': "Cover",
        '~': "Independence",
        '#': "Mutual Exclusion",
    }

    # MacCartney's join table for composing relations (B o A)
    # This is a partial table containing only the entries needed for this problem.
    JOIN_TABLE = {
        'sq': {  # Relation A
            '#': '#'  # Relation B: join(sq, #) -> #
        },
        '#': {   # Relation A
            '|': 'sq'  # Relation B: join(#, |) -> sq
        }
    }

    # Step 1: Deconstruct the problem into a path: Premise -> Intermediate -> Hypothesis
    premise = "Mark is singing a pop song by Taylor Swift"
    intermediate = "Mark is singing a song by Michael Jackson"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    print("Step 1: Deconstructing the inference into a path P -> I -> H.")
    print(f"  P (Premise):      '{premise}'")
    print(f"  I (Intermediate): '{intermediate}'")
    print(f"  H (Hypothesis):   '{hypothesis}'\n")

    # Step 2: Determine the relation from Premise (P) to Intermediate (I)
    print("Step 2: Determining rel(P, I) by composing lexical relations.")
    
    # Lexical relation 1: "pop song" vs "song"
    rel_lex1_sym = 'sq'
    # Lexical relation 2: "Taylor Swift" vs "Michael Jackson"
    rel_lex2_sym = '#'
    
    print(f"  Relation('pop song', 'song') is '{RELATIONS[rel_lex1_sym]}' ({rel_lex1_sym}).")
    print(f"  Relation('Taylor Swift', 'Michael Jackson') is '{RELATIONS[rel_lex2_sym]}' ({rel_lex2_sym}).")

    # Compose the lexical relations to find rel(P, I)
    rel_P_I_sym = JOIN_TABLE[rel_lex1_sym][rel_lex2_sym]
    print(f"  Composition: join({rel_lex1_sym}, {rel_lex2_sym}) = {rel_P_I_sym}")
    print(f"  Therefore, rel(P, I) is '{RELATIONS[rel_P_I_sym]}' ({rel_P_I_sym}).\n")

    # Step 3: Determine the relation from Intermediate (I) to Hypothesis (H)
    print("Step 3: Determining rel(I, H).")
    # H is the negation of I
    rel_I_H_sym = '|'
    print(f"  H is the negation of I, so rel(I, H) is '{RELATIONS[rel_I_H_sym]}' ({rel_I_H_sym}).\n")

    # Step 4: Compute the final relation from P to H
    print("Step 4: Computing the final relation rel(P, H) by composition.")
    # Final composition: join(rel(P, I), rel(I, H))
    final_rel_sym = JOIN_TABLE[rel_P_I_sym][rel_I_H_sym]
    final_rel_name = RELATIONS[final_rel_sym]

    print(f"  Final Composition: join(rel(P, I), rel(I, H)) => join({rel_P_I_sym}, {rel_I_H_sym}) = {final_rel_sym}")
    print("----------------------------------------------------------------")
    print(f"The name of the final projected natural logic operator is: {final_rel_name}")
    print("----------------------------------------------------------------")

if __name__ == '__main__':
    solve_maccartney_inference()