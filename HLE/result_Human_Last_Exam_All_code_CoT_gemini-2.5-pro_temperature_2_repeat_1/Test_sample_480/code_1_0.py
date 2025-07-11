def solve_mccartney_inference():
    """
    Solves for the final natural logic operator for the given P/H pair
    using MacCartney's compositional proof framework.
    """
    
    # The 7 natural logic relations
    RELATIONS = {
        'eq': {'name': 'Identity', 'symbol': '='},
        'fwd_entail': {'name': 'Forward Entailment', 'symbol': '⊂'},
        'rev_entail': {'name': 'Reverse Entailment', 'symbol': '⊃'},
        'negation': {'name': 'Negation', 'symbol': '^'},
        'alternation': {'name': 'Alternation', 'symbol': '|'},
        'cover': {'name': 'Cover', 'symbol': '∪'},
        'independence': {'name': 'Independence', 'symbol': '~'}
    }
    
    # MacCartney's composition table (Join Table `R_old ⨢ R_edit`)
    # Rows are the previous relation, columns are the relation of the current edit
    # Let's use 'fe' for fwd_entail, 're' for rev_entail for brevity
    composition_table = {
        #             eq,     fe,       re,           neg,         alt,          cov,          ind
        'eq':        ['eq', 'fwd_entail', 'rev_entail', 'negation', 'alternation', 'cover', 'independence'],
        'fwd_entail':['fwd_entail','fwd_entail','independence','negation', 'alternation', 'cover', 'independence'],
        'rev_entail':['rev_entail','independence','rev_entail','cover',    'alternation', 'cover', 'independence'],
        'negation':  ['negation','alternation','cover',     'eq',       'fwd_entail', 'rev_entail', 'independence'],
        'alternation':['alternation','fwd_entail','rev_entail', 'fwd_entail', 'eq',       'rev_entail', 'independence'],
        'cover':     ['cover',     'rev_entail','independence','rev_entail', 'fwd_entail', 'eq', 'independence'],
        'independence':['independence','independence','independence','independence','independence','independence','independence']
    }
    col_map = {name: i for i, name in enumerate(['eq', 'fwd_entail', 'rev_entail', 'negation', 'alternation', 'cover', 'independence'])}

    # Projection of a lexical relation (row) through a context (col)
    # Contexts: = (upward monotone), ^ (downward monotone)
    projection_table = {
        #             '=' (ID)    '^' (NEG)
        'fwd_entail': ['fwd_entail', 'rev_entail'],
        'rev_entail': ['rev_entail', 'fwd_entail'],
        'negation':   ['negation', 'eq'],
        'alternation':['alternation', 'cover'],
        'cover':      ['cover', 'alternation'],
        'independence': ['independence', 'independence']
    }

    # -- Simulation of the Proof --
    
    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    print(f"Premise: \"{premise}\"")
    print(f"Hypothesis: \"{hypothesis}\"")
    print("-" * 30)

    current_relation = 'eq'
    print(f"Step 0: Initial Relation. The premise entails itself.")
    print(f"         P => P. Relation is {RELATIONS[current_relation]['name']} ({RELATIONS[current_relation]['symbol']})")
    print("-" * 30)
    
    # Edit 1: Insertion of "not"
    lexical_relation_1 = 'negation'
    context_1 = 'eq' # The initial context is positive / identity
    projected_edit_1 = lexical_relation_1
    
    print(f"Step 1: Insert 'not'.")
    print(f"         Lexical Relation: {RELATIONS[lexical_relation_1]['name']} ({RELATIONS[lexical_relation_1]['symbol']})")
    old_relation = current_relation
    current_relation = composition_table[old_relation][col_map[projected_edit_1]]
    print(f"         Composition: {RELATIONS[old_relation]['name']} ⨢ {RELATIONS[projected_edit_1]['name']} = {RELATIONS[current_relation]['name']} ({RELATIONS[current_relation]['symbol']})")
    print(f"         New Projected Relation is {RELATIONS[current_relation]['name']} ({RELATIONS[current_relation]['symbol']})")
    print("-" * 30)

    # Edit 2: Deletion of "pop"
    lexical_relation_2 = 'fwd_entail'
    context_2 = 'negation' # The context is now negative because of 'not'
    projected_edit_2 = projection_table[lexical_relation_2][1]
    
    print(f"Step 2: Delete 'pop'.")
    print(f"         Lexical Relation ('pop song' to 'song'): {RELATIONS[lexical_relation_2]['name']} ({RELATIONS[lexical_relation_2]['symbol']})")
    print(f"         Context is negative, so the relation is projected from ⊂ to ⊃.")
    print(f"         Projected Edit Relation: {RELATIONS[projected_edit_2]['name']} ({RELATIONS[projected_edit_2]['symbol']})")
    old_relation = current_relation
    current_relation = composition_table[old_relation][col_map[projected_edit_2]]
    print(f"         Composition: {RELATIONS[old_relation]['name']} ⨢ {RELATIONS[projected_edit_2]['name']} = {RELATIONS[current_relation]['name']} ({RELATIONS[current_relation]['symbol']})")
    print(f"         New Projected Relation is {RELATIONS[current_relation]['name']} ({RELATIONS[current_relation]['symbol']})")
    print("-" * 30)
    
    # Edit 3: Substitution of "Taylor Swift" -> "Michael Jackson"
    lexical_relation_3 = 'alternation'
    context_3 = 'negation' # Context remains negative
    projected_edit_3 = projection_table[lexical_relation_3][1]
    
    print(f"Step 3: Substitute 'Taylor Swift' with 'Michael Jackson'.")
    print(f"         Lexical Relation: {RELATIONS[lexical_relation_3]['name']} ({RELATIONS[lexical_relation_3]['symbol']})")
    print(f"         Context is negative, so the relation is projected from | to ∪.")
    print(f"         Projected Edit Relation: {RELATIONS[projected_edit_3]['name']} ({RELATIONS[projected_edit_3]['symbol']})")
    old_relation = current_relation
    current_relation = composition_table[old_relation][col_map[projected_edit_3]]
    print(f"         Composition: {RELATIONS[old_relation]['name']} ⨢ {RELATIONS[projected_edit_3]['name']} = {RELATIONS[current_relation]['name']} ({RELATIONS[current_relation]['symbol']})")
    print(f"         New Projected Relation is {RELATIONS[current_relation]['name']} ({RELATIONS[current_relation]['symbol']})")
    print("-" * 30)

    final_operator_name = RELATIONS[current_relation]['name']
    print(f"The final projected natural logic operator is: {final_operator_name}")

solve_mccartney_inference()