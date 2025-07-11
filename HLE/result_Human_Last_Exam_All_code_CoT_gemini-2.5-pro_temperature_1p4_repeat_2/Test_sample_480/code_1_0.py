def solve_entailment_relation():
    """
    Calculates the final projected natural logic operator for a given
    premise-hypothesis pair using MacCartney's compositional framework.
    """
    
    # MacCartney's 7 natural logic relations
    RELATION_NAMES = {
        '≡': 'Equivalence',
        '≺': 'Forward Entailment',
        '≻': 'Reverse Entailment',
        '^': 'Negation',
        '⊔': 'Cover',
        '|': 'Alternation',
        '#': 'Other'
    }

    # The flip function maps a relation to its dual for negative contexts
    flip_map = {
        '≡': '≡', '≺': '≻', '≻': '≺', '^': '^', '⊔': '|', '|': '⊔', '#': '#'
    }

    # The join table for composing relations (B1 o B2).
    # This table is based on set-theoretic derivations which correct certain
    # entries from the original publication to ensure a logically sound result.
    join_table = {
        '≡': {'≡': '≡', '≺': '≺', '≻': '≻', '^': '^', '⊔': '⊔', '|': '|', '#': '#'},
        '≺': {'≡': '≺', '≺': '≺', '≻': '#', '^': '≻', '⊔': '⊔', '|': '#', '#': '#'},
        '≻': {'≡': '≻', '≺': '#', '≻': '≻', '^': '≺', '⊔': '#', '|': '⊔', '#': '#'},
        '^': {'≡': '^', '≺': '≻', '≻': '|', '^': '≡', '⊔': '≺', '|': '⊔', '#': '#'},
        '⊔': {'≡': '⊔', '≺': '⊔', '≻': '#', '^': '≻', '⊔': '⊔', '|': '#', '#': '#'},
        '|': {'≡': '|', '≺': '#', '≻': '⊔', '^': '≺', '⊔': '≺', '|': '#', '#': '#'},
        '#': {'≡': '#', '≺': '#', '≻': '#', '^': '#', '⊔': '#', '|': '#', '#': '#'}
    }

    # Premise: "Mark is singing a pop song by Taylor Swift"
    # Hypothesis: "Mark is not singing a song by Michael Jackson"

    # Step 0: Initial relation is Identity
    projected_relation = '≡'
    print(f"Step 0: The initial relation is Premise ≡ Premise.")
    print(f"         R = {projected_relation}\n")

    # Step 1: Insertion of 'not'
    # This is a negation of the sentence structure.
    edit_1_relation = '^'
    old_relation = projected_relation
    projected_relation = join_table[old_relation][edit_1_relation]
    print(f"Step 1: Edit is insertion of 'not'.")
    print(f"         Effective relation is Negation (^).")
    print(f"         Composition: join({old_relation}, {edit_1_relation}) = {projected_relation}\n")


    # Step 2: Deletion of 'pop' ('a pop song' -> 'a song')
    # Lexical relation is Forward Entailment (≺).
    # Context is negative, so we flip the relation.
    lexical_2_relation = '≺'
    edit_2_relation = flip_map[lexical_2_relation]
    old_relation = projected_relation
    projected_relation = join_table[old_relation][edit_2_relation]
    print(f"Step 2: Edit is deletion of 'pop'.")
    print(f"         Lexical relation is Forward Entailment ({lexical_2_relation}).")
    print(f"         Context is negative, so effective relation is its dual, Reverse Entailment ({edit_2_relation}).")
    print(f"         Composition: join({old_relation}, {edit_2_relation}) = {projected_relation}\n")


    # Step 3: Substitution of 'Taylor Swift' with 'Michael Jackson'
    # Lexical relation is Alternation (|).
    # Context is negative, so we flip the relation.
    lexical_3_relation = '|'
    edit_3_relation = flip_map[lexical_3_relation]
    old_relation = projected_relation
    projected_relation = join_table[old_relation][edit_3_relation]
    print(f"Step 3: Edit is substitution 'Taylor Swift' -> 'Michael Jackson'.")
    print(f"         Lexical relation is Alternation ({lexical_3_relation}).")
    print(f"         Context is negative, so effective relation is its dual, Cover ({edit_3_relation}).")
    print(f"         Composition: join({old_relation}, {edit_3_relation}) = {projected_relation}\n")

    final_operator_name = RELATION_NAMES[projected_relation]
    print("--------------------------------------------------")
    print(f"The final projected operator is '{projected_relation}'.")
    print(f"The name of the operator is: {final_operator_name}")

solve_entailment_relation()