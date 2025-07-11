def solve_nli_macartney():
    """
    Solves the entailment problem using MacCartney's natural logic framework.
    """
    # MacCartney's 7 relations:
    # = : equivalence
    # [ : forward entailment (subset)
    # ] : reverse entailment (superset)
    # ^ : negation (complement)
    # | : cover (exhaustion)
    # v : alternation (disjointness)
    # # : independence

    relations = ['=', '[', ']', '^', '|', 'v', '#']
    relation_names = {
        '=': 'Equivalence',
        '[': 'Forward Entailment',
        ']': 'Reverse Entailment',
        '^': 'Negation',
        '|': 'Cover',
        'v': 'Alternation',
        '#': 'Independence',
    }

    # MacCartney's projection join table B(atomic, context)
    # Some indeterminate cells are filled with common outcomes (#, |, etc.)
    B_table_str = """
    = [ ] ^ | v #
    [ [ # v | v #
    ] # ] | v | #
    ^ ] [ = [ ] #
    | ] v [ | # #
    v # # ] # # #
    # # # # # # #
    """
    
    # Parse the table into a dictionary for easy lookup
    b_table = {}
    rows = B_table_str.strip().split('\n')
    row_headers = ['=', '[', ']', '^', '|', 'v', '#']
    col_headers = ['=', '[', ']', '^', '|', 'v', '#']
    
    for i, row_str in enumerate(rows):
        row_header = row_headers[i]
        b_table[row_header] = {}
        parts = row_str.strip().split()
        for j, cell in enumerate(parts):
            col_header = col_headers[j]
            b_table[row_header][col_header] = cell


    # Premise: "Mark is singing a pop song by Taylor Swift"
    # Hypothesis: "Mark is not singing a song by Michael Jackson"

    # Edits executed from left-to-right on the hypothesis structure
    edits = [
        {"desc": "Insert 'not'", "rel": "^"},
        {"desc": "Delete 'pop' (pop song -> song)", "rel": "]"},
        {"desc": "Substitute 'Taylor Swift' with 'Michael Jackson'", "rel": "v"},
    ]

    # Start with the identity relation
    current_relation = '='
    print("Starting entailment proof...")
    print(f"Initial Relation: {current_relation} ({relation_names[current_relation]})")
    print("-" * 30)

    # Compose the relations sequentially
    for i, edit in enumerate(edits):
        atomic_rel = edit["rel"]
        print(f"Step {i+1}: Applying Edit '{edit['desc']}'")
        print(f"  Atomic Relation: {atomic_rel} ({relation_names[atomic_rel]})")
        print(f"  Current Context Relation: {current_relation} ({relation_names[current_relation]})")
        
        # Look up the new relation in the projection table B(atomic, context)
        # B(atomic, current)
        new_relation = b_table[atomic_rel][current_relation]
        
        print(f"  Composition: B({atomic_rel}, {current_relation}) = {new_relation}")
        current_relation = new_relation
        print(f"  New Projected Relation: {current_relation} ({relation_names[current_relation]})")
        print("-" * 30)

    final_relation_name = relation_names[current_relation]
    print(f"Final Projected Natural Logic Operator: {final_relation_name}")

solve_nli_macartney()