import collections

def solve_entailment():
    """
    Calculates the final projected natural logic operator for the given inference
    using MacCartney's compositional framework.
    """
    # Define the 7 semantic relations and their symbols for convenience
    relations = {
        'eq': ('=', 'Identity'),
        'f_ent': ('⊏', 'Forward Entailment'),
        'r_ent': ('⊐', 'Reverse Entailment'),
        'neg': ('^', 'Negation'),
        'alt': ('|', 'Alternation'),
        'cov': ('#', 'Cover'),
        'ind': ('~', 'Independence')
    }

    # MacCartney's composition table
    # comp_table[rel1][rel2] gives the composition of rel1 ; rel2
    comp_table = {
        'eq':    {'eq': 'eq', 'f_ent': 'f_ent', 'r_ent': 'r_ent', 'neg': 'neg', 'alt': 'alt', 'cov': 'cov', 'ind': 'ind'},
        'f_ent': {'eq': 'f_ent', 'f_ent': 'f_ent', 'r_ent': 'ind', 'neg': 'alt', 'alt': 'alt', 'cov': 'ind', 'ind': 'ind'},
        'r_ent': {'eq': 'r_ent', 'f_ent': 'ind', 'r_ent': 'r_ent', 'neg': 'cov', 'alt': 'ind', 'cov': 'cov', 'ind': 'ind'},
        'neg':   {'eq': 'neg', 'f_ent': 'cov', 'r_ent': 'alt', 'neg': 'eq', 'alt': 'r_ent', 'cov': 'f_ent', 'ind': 'ind'},
        'alt':   {'eq': 'alt', 'f_ent': 'ind', 'r_ent': 'alt', 'neg': 'f_ent', 'alt': 'ind', 'cov': 'f_ent', 'ind': 'ind'},
        'cov':   {'eq': 'cov', 'f_ent': 'cov', 'r_ent': 'ind', 'neg': 'r_ent', 'alt': 'r_ent', 'cov': 'ind', 'ind': 'ind'},
        'ind':   {'eq': 'ind', 'f_ent': 'ind', 'r_ent': 'ind', 'neg': 'ind', 'alt': 'ind', 'cov': 'ind', 'ind': 'ind'}
    }

    # Premise: "Mark is singing a pop song by Taylor Swift"
    # Hypothesis: "Mark is not singing a song by Michael Jackson"
    # Edits are processed left-to-right based on the hypothesis structure.
    
    edits = [
        {
            "description": "Insertion of 'not'",
            "lexical_relation": "neg",
            "context": "None (atomic edit)",
            "projected_relation": "neg"
        },
        {
            "description": "Deletion of 'pop' ('pop song' -> 'song')",
            "lexical_relation": "f_ent", # 'pop song' ⊏ 'song'
            "context": "Downward-monotone (inside 'not')",
            "projected_relation": "r_ent" # ⊏ projects to ⊐ in ↓ context
        },
        {
            "description": "Substitution 'Taylor Swift' -> 'Michael Jackson'",
            "lexical_relation": "alt", # 'Taylor Swift' | 'Michael Jackson'
            "context": "Downward-monotone (inside 'not')",
            "projected_relation": "cov" # | projects to # in ↓ context
        }
    ]

    current_relation = 'eq' # Start with the Identity relation

    print("Premise: 'Mark is singing a pop song by Taylor Swift'")
    print("Hypothesis: 'Mark is not singing a song by Michael Jackson'\n")
    print("Beginning entailment calculation by composing relations from each edit...")
    print("------------------------------------------------------------------")

    for i, edit in enumerate(edits):
        previous_relation = current_relation
        projected_rel = edit["projected_relation"]
        current_relation = comp_table[previous_relation][projected_rel]

        prev_sym, prev_name = relations[previous_relation]
        proj_sym, proj_name = relations[projected_rel]
        curr_sym, curr_name = relations[current_relation]

        print(f"Step {i+1}: {edit['description']}")
        print(f"  - Lexical Relation: {relations[edit['lexical_relation']][1]} ({relations[edit['lexical_relation']][0]})")
        print(f"  - Context: {edit['context']}")
        print(f"  - Projected Relation: {proj_name} ({proj_sym})")
        print(f"  - Composition: {prev_name} ; {proj_name} => {curr_name}")
        print(f"  - Symbolic Equation: {prev_sym} ; {proj_sym} = {curr_sym}")
        print("------------------------------------------------------------------")

    final_symbol, final_name = relations[current_relation]
    print(f"\nFinal Result:")
    print(f"The final projected natural logic operator is '{final_name}', represented by the symbol '{final_symbol}'.")

solve_entailment()
<<<Forward Entailment>>>