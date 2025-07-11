import sys
from io import StringIO

def solve():
    """
    Calculates the final projected natural logic operator for the given inference
    using MacCartney's compositional semantics.
    """
    RELATIONS = {
        '≡': 'equivalence',
        '⊂': 'forward entailment',
        '⊃': 'reverse entailment',
        '^': 'negation',
        '|': 'alternation',
        'v': 'cover',
        '#': 'independence',
    }

    # Composition table based on MacCartney (2009), Table 3.5
    # Rows are the current relation, columns are the projected relation of the new edit.
    COMPOSITION_TABLE = {
        #      ≡    ⊂    ⊃    ^    |    v    #
        '≡': {'≡':'≡', '⊂':'⊂', '⊃':'⊃', '^':'^', '|':'|', 'v':'v', '#':'#'},
        '⊂': {'≡':'⊂', '⊂':'⊂', '⊃':'#', '^':'^', '|':'|', 'v':'#', '#':'#'},
        '⊃': {'≡':'⊃', '⊂':'#', '⊃':'⊃', '^':'#', '|':'|', 'v':'v', '#':'v'},
        '^': {'≡':'^', '⊂':'v', '⊃':'|', '^':'≡', '|':'⊃', 'v':'⊂', '#':'#'},
        '|': {'≡':'|', '⊂':'#', '⊃':'v', '^':'v', '|':'|', 'v':'⊂', '#':'⊂'},
        'v': {'≡':'v', '⊂':'⊃', '⊃':'#', '^':'⊂', '|':'⊂', 'v':'v', '#':'#'},
        '#': {'≡':'#', '⊂':'#', '⊃':'#', '^':'#', '|':'|', 'v':'v', '#':'#'},
    }
    
    # Projection table based on MacCartney (2009), Table 3.4
    PROJECTION_TABLE = {
        #             ≡    ⊂    ⊃    ^    |    v    #
        'UPWARD':   {'≡':'≡', '⊂':'⊂', '⊃':'⊃', '^':'^', '|':'|', 'v':'v', '#':'#'},
        'DOWNWARD': {'≡':'≡', '⊂':'⊃', '⊃':'⊂', '^':'^', '|':'v', 'v':'|', '#':'#'},
    }

    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    # Define the edits required to transform the premise into the hypothesis.
    edits = [
        {
            "description": "Insert 'not'",
            "lexical_relation": '|',
            "context": "UPWARD" # The word is inserted into a positive context
        },
        {
            "description": "Delete 'pop'",
            "lexical_relation": '⊂', # 'pop song' is a subtype of 'song'
            "context": "DOWNWARD" # Occurs under the scope of 'not'
        },
        {
            "description": "Substitute 'Taylor Swift' with 'Michael Jackson'",
            "lexical_relation": 'v', # 'TS' and 'MJ' are disjoint sets (A ⊂ ¬B)
            "context": "DOWNWARD" # Occurs under the scope of 'not'
        }
    ]

    print("Deriving the entailment relation for:")
    print(f'Premise:  "{premise}"')
    print(f'Hypothesis: "{hypothesis}"\n')

    # Start with equivalence, since P entails P
    current_relation = '≡'
    print(f"Starting with initial relation: {current_relation} ({RELATIONS[current_relation]})\n")

    for i, edit in enumerate(edits):
        lex_rel = edit['lexical_relation']
        context = edit['context']

        # 1. Project the lexical relation through its semantic context
        projected_relation = PROJECTION_TABLE[context][lex_rel]

        # 2. Compose the running relation with the projected relation
        new_relation = COMPOSITION_TABLE[current_relation][projected_relation]

        # Print the derivation step
        print(f"Step {i+1}: {edit['description']}")
        print(f"  - Lexical relation: {lex_rel} ({RELATIONS[lex_rel]})")
        print(f"  - Semantic context: {context.lower()}")
        print(f"  - Projected relation: {projected_relation} ({RELATIONS[projected_relation]})")
        # Final output line must show the final equation
        print(f"  - Composition: {current_relation} ○ {projected_relation} = {new_relation}")
        print(f"  - New relation: {new_relation} ({RELATIONS[new_relation]})\n")
        
        current_relation = new_relation

    final_relation_symbol = current_relation
    final_relation_name = RELATIONS[final_relation_symbol]

    print("---" * 10)
    print(f"Final Projected Relation: {final_relation_symbol}")
    print(f"Name of Operator: {final_relation_name}")

solve()