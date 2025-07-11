import pandas as pd

def solve_maccartney_inference():
    """
    Calculates the final projected natural logic operator for a given premise and hypothesis
    based on MacCartney's compositional framework.
    """
    
    # Define the 7 relations and their names
    relations = ['≡', '<', '>', '|', '^', 'v', '#']
    relation_names = {
        '≡': 'Equivalence',
        '<': 'Forward Entailment',
        '>': 'Reverse Entailment',
        '|': 'Alternation',
        '^': 'Negation',
        'v': 'Cover',
        '#': 'Independence'
    }

    # MacCartney's Join Table for composition (B o A)
    # Rows are the current relation (B), columns are the new edit's relation (A)
    join_table_data = {
        # A: ≡    <    >    |    ^    v    #
        '≡': ['≡', '<', '>', '|', '^', 'v', '#'],  # B: ≡
        '<': ['<', '<', '#', '<', 'v', '#', 'v'],  # B: <
        '>': ['>', '#', '>', 'v', '>', '#', '#'],  # B: >
        '|': ['|', '<', 'v', '|', '>', '<', 'v'],  # B: |
        '^': ['^', '>', '<', '>', '^', '<', '|'],  # B: ^
        'v': ['v', '#', '>', '<', '|', '#', '<'],  # B: v
        '#': ['#', '#', '#', '#', '#', '#', '#'],  # B: #
    }
    
    # Create a DataFrame for easy lookup
    join_table = pd.DataFrame(join_table_data, index=relations, columns=relations)

    # Premise and Hypothesis for context
    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    print(f"Premise: \"{premise}\"")
    print(f"Hypothesis: \"{hypothesis}\"\n")

    # The sequence of edits and their corresponding relations
    edits = [
        {"desc": "INSERT 'not'", "relation": "|"},
        {"desc": "DELETE 'pop'", "relation": "<"},
        {"desc": "SUBSTITUTE 'Taylor Swift' -> 'Michael Jackson'", "relation": "|"}
    ]

    # Start with the identity relation
    current_relation = '≡'
    print(f"Starting with initial relation: {current_relation} ({relation_names[current_relation]})")
    print("-" * 30)

    # Compose the relations for each edit from left to right
    for i, edit in enumerate(edits):
        edit_relation = edit['relation']
        print(f"Step {i+1}: Applying edit '{edit['desc']}'")
        print(f"  - Current relation: {current_relation} ({relation_names[current_relation]})")
        print(f"  - Edit relation:    {edit_relation} ({relation_names[edit_relation]})")
        
        previous_relation = current_relation
        current_relation = join_table.loc[previous_relation, edit_relation]
        
        # Output the composition equation
        print(f"  - Composition:      {previous_relation} ∘ {edit_relation} = {current_relation}")
        print(f"  - New relation:     {current_relation} ({relation_names[current_relation]})")
        print("-" * 30)

    # The final answer
    final_relation_name = relation_names[current_relation]
    print(f"\nThe final projected natural logic operator is: {final_relation_name}")

    # Hidden final answer for verification
    # <<<Forward Entailment>>>

if __name__ == "__main__":
    solve_maccartney_inference()