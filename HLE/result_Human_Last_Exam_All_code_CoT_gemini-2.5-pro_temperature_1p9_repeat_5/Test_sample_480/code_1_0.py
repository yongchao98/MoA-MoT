def solve_mccartney_inference():
    """
    Solves for the final projected natural logic operator for the given P/H pair
    by implementing MacCartney's compositional semantics.
    """
    
    # The 7 semantic relations and their names
    RELATION_NAMES = {
        '=': 'Equivalence',
        'sq': 'Forward Entailment',
        'rev': 'Reverse Entailment',
        '^': 'Negation (Contradiction)',
        'v': 'Alternation (Exhaustive Opposition)',
        '|': 'Cover (Non-Contradiction)',
        '#': 'Independence',
    }

    # The join table for composing relations (r_c = row, r_a = column)
    # Source: Bowman, Angeli, Potts, MacCartney (2015)
    JOIN_TABLE = {
        '=':   {'=': '=',  'sq': 'sq',  'rev': 'rev', '^': '^', 'v': 'v', '|': '|', '#': '#'},
        'sq':  {'=': 'sq', 'sq': 'sq',  'rev': '#',   '^': 'sq', 'v': 'v', '|': '#', '#': '#'},
        'rev': {'=': 'rev','sq': '#',   'rev': 'rev', '^': '#', 'v': 'v', '|': 'rev', '#': '#'},
        '^':   {'=': '^',  'sq': '#',   'rev': '^',   '^': '#', 'v': '|', '|': 'sq', '#': '#'},
        'v':   {'=': 'v',  'sq': 'v',   'rev': '#',   '^': 'v', 'v': '=', '|': 'sq', '#': '#'},
        '|':   {'=': '|',  'sq': 'sq',  'rev': '|',   '^': 'sq', 'v': '=', '|': '|', '#': '#'},
        '#':   {'=': '#',  'sq': '#',   'rev': '#',   '^': '#', 'v': '#', '|': '#', '#': '#'},
    }

    # Relations that create downward-monotone (negating) contexts
    DOWNWARD_MONOTONE_CONTEXTS = {'rev', '^', 'v', '|'}

    # How relations flip in downward-monotone contexts
    FLIP_RULES = {
        'sq': 'rev',
        'rev': 'sq',
        '^': '|',
        '|': '^',
    }

    # The sequence of edits and their atomic relations, processed left-to-right
    edits = [
        {"desc": "Insert 'not'", "rel": "v"},
        {"desc": "Delete 'pop'", "rel": "rev"},
        {"desc": "Substitute 'Taylor Swift' -> 'Michael Jackson'", "rel": "^"},
    ]

    # --- Composition Logic ---
    projected_relation = '='
    equation_steps = [projected_relation]
    
    print("Premise: 'Mark is singing a pop song by Taylor Swift'")
    print("Hypothesis: 'Mark is not singing a song by Michael Jackson'")
    print("-" * 40)
    print(f"Starting with initial relation: {projected_relation} ({RELATION_NAMES[projected_relation]})")

    for i, edit in enumerate(edits):
        print("-" * 40)
        print(f"Step {i+1}: Applying edit '{edit['desc']}'")
        
        atomic_relation = edit['rel']
        print(f"  Atomic edit relation: {atomic_relation} ({RELATION_NAMES[atomic_relation]})")

        context_relation = projected_relation
        print(f"  Current context relation is {context_relation} ({RELATION_NAMES[context_relation]})")

        is_downward = context_relation in DOWNWARD_MONOTONE_CONTEXTS
        if is_downward:
            print(f"  Context is Downward Monotone.")
            final_edit_relation = FLIP_RULES.get(atomic_relation, atomic_relation)
            if final_edit_relation != atomic_relation:
                print(f"  Flipping atomic relation '{atomic_relation}' -> '{final_edit_relation}'")
                equation_steps.append(f"flip({atomic_relation})")
            else:
                equation_steps.append(atomic_relation)
        else:
            print("  Context is Upward Monotone. Relation does not flip.")
            final_edit_relation = atomic_relation
            equation_steps.append(atomic_relation)
        
        # Compose the new relation
        new_projected_relation = JOIN_TABLE[context_relation][final_edit_relation]
        print(f"  Composition: join({context_relation}, {final_edit_relation}) = {new_projected_relation}")
        projected_relation = new_projected_relation
    
    print("-" * 40)
    final_name = RELATION_NAMES[projected_relation]
    print(f"Final projected relation is '{projected_relation}', which stands for: {final_name}")
    
    # Constructing the final equation string
    # e.g., join(join(join(=, v), flip(rev)), flip(^))
    equation = equation_steps.pop(0) # start with '='
    while equation_steps:
      equation = f"join({equation}, {equation_steps.pop(0)})"
    
    print(f"\nFinal Equation: {projected_relation} = {equation}")


solve_mccartney_inference()