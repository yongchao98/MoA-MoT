import collections

def solve_mccartney_inference():
    """
    Solves the MacCartney natural logic inference problem step-by-step.
    """

    # 1. Define the premise, hypothesis, and transformation steps.
    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    # Define the sequence of edits from left to right.
    edits = [
        {
            "description": "Delete 'pop' from 'a pop song'",
            "intermediate_result": "Mark is singing a song by Taylor Swift",
            "relation_symbol": "<",
            "relation_name": "Forward Entailment"
        },
        {
            "description": "Insert 'not' to negate 'is singing'",
            "intermediate_result": "Mark is not singing a song by Taylor Swift",
            "relation_symbol": "^",
            "relation_name": "Negation"
        },
        {
            "description": "Substitute 'Taylor Swift' with 'Michael Jackson'",
            "intermediate_result": "Mark is not singing a song by Michael Jackson",
            "relation_symbol": "v",
            "relation_name": "Cover"
        }
    ]

    # 2. Define the relational symbols and MacCartney's composition table.
    relations = ['=', '<', '>', '^', '|', 'v', '#']
    rel_map = {name: i for i, name in enumerate(relations)}
    rel_name_map = {
        '=': 'Equivalence',
        '<': 'Forward Entailment',
        '>': 'Reverse Entailment',
        '^': 'Negation',
        '|': 'Alternation',
        'v': 'Cover',
        '#': 'Independence'
    }

    # Composition table: `table[row][col]` is the composition of `row ; col`
    composition_table_str = """
      = < > ^ | v #
      < < # ^ | < # #
      > # > v # > #
      ^ ^ v = v ^ #
      | v < > # # # #
      v # | < | # v
      # # # # # # # #
    """
    rows = composition_table_str.strip().split('\n')
    table = [row.split() for row in rows]
    
    # 3. Perform the compositional inference.
    print(f"Premise: \"{premise}\"\n")
    print("Starting inference process...\n")

    current_relation_symbol = "="
    current_sentence = premise

    # Initial state
    print(f"Step 0: Initial State")
    print(f"  - Sentence: \"{current_sentence}\"")
    print(f"  - Current projected relation: {current_relation_symbol} ({rel_name_map[current_relation_symbol]})\n")
    
    for i, edit in enumerate(edits):
        print(f"Step {i + 1}: {edit['description']}")
        
        edit_relation_symbol = edit['relation_symbol']
        print(f"  - Atomic edit relation: {edit_relation_symbol} ({rel_name_map[edit_relation_symbol]})")
        
        # Look up the composition in the table
        row_idx = rel_map[current_relation_symbol]
        col_idx = rel_map[edit_relation_symbol]
        new_relation_symbol = table[row_idx][col_idx]
        
        print(f"  - Composition: {rel_name_map[current_relation_symbol]} ({current_relation_symbol}) ; {rel_name_map[edit_relation_symbol]} ({edit_relation_symbol}) => {rel_name_map[new_relation_symbol]} ({new_relation_symbol})")
        
        current_relation_symbol = new_relation_symbol
        current_sentence = edit['intermediate_result']
        
        print(f"  - Intermediate Sentence: \"{current_sentence}\"")
        print(f"  - New projected relation: {current_relation_symbol} ({rel_name_map[current_relation_symbol]})\n")

    print("--------------------------------------------------")
    print(f"Final Hypothesis: \"{hypothesis}\"")
    print(f"The final projected natural logic operator is: {rel_name_map[current_relation_symbol]}")
    print("--------------------------------------------------")


solve_mccartney_inference()