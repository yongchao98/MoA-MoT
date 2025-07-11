import sys

def solve_natural_logic_inference():
    """
    Calculates the final projected natural logic operator for a given
    premise-hypothesis pair based on MacCartney's framework.
    """
    # Define the names for the 7 semantic relations (+ 'unknown')
    relation_names = {
        '=': 'Equivalence',
        '<': 'Forward Entailment',
        '>': 'Reverse Entailment',
        '^': 'Negation',
        '|': 'Alternation',
        'v': 'Cover', # 'v' is also called Union or Join
        '#': 'Independence', # '#' is also called Cover in some papers. We follow the 2009 paper where '#' is Independence/Exhaustive. Wait, no, lets check the source paper.
                           # Let's clarify based on MacCartney's 2009 paper "Natural Language Inference":
                           # | : alternation (mutually exclusive), e.g., cat | dog
                           # v : cover (union, exhaustive), e.g., animal v mammal
                           # # : independence, e.g., cat # hungry.
                           # The relationship between "song by Taylor Swift" and "song by Michael Jackson"
                           # is mutual exclusion, so | (Alternation) is the best fit.
                           # Let's re-calculate with '|'
                           # Step 1: = o | = |
                           # Step 2: | o ^ = <
                           # The result is the same! But using the correct term is better.
        '#': 'Independence',
        'u': 'Unknown',
    }
    relation_names['|'] = 'Alternation' # Overriding for clarity in this problem.

    # MacCartney's composition table (B_row ◦ A_col)
    # Source: MacCartney, "Natural Language Inference", 2009.
    composition_table = {
        '=': {'=': '=', '<': '<', '>': '>', '^': '^', '|': '|', 'v': 'v', '#': '#'},
        '<': {'=': '<', '<': '<', '>': 'u', '^': '|', '|': 'u', 'v': '<', '#': 'u'},
        '>': {'=': '>', '<': 'v', '>': '>', '^': '#', '|': '>', 'v': 'u', '#': 'u'},
        '^': {'=': '^', '<': '>', '>': '<', '^': '=', '|': 'v', 'v': '|', '#': '#'},
        '|': {'=': '|', '<': 'u', '>': 'v', '^': '<', '|': 'u', 'v': '<', '#': 'u'},
        'v': {'=': 'v', '<': 'v', '>': 'u', '^': '>', '|': '>', 'v': 'u', '#': 'u'},
        '#': {'=': '#', '<': 'u', '>': '>', '^': '<', '|': 'u', 'v': '>', '#': 'u'},
    }
    
    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"
    
    print(f"Premise (P): \"{premise}\"")
    print(f"Hypothesis (H): \"{hypothesis}\"")
    print("-" * 25)
    print("Executing edits from left to right to transform P into H.")
    
    # --- Initial State ---
    current_relation = '='
    print(f"\nStep 0: Initial State")
    print(f"  The relation of P to itself is {current_relation} ({relation_names[current_relation]})")

    # --- Step 1: Substitution ---
    intermediate_text = "Mark is singing a song by Michael Jackson"
    # The relationship between "singing a song by TS" and "singing a song by MJ" is
    # that they are mutually exclusive. This is Alternation.
    edit1_relation = '|' 
    print(f"\nStep 1: Edit Premise")
    print(f"  Action: Substitute 'a pop song by Taylor Swift' with 'a song by Michael Jackson'")
    print(f"  Intermediate Text: \"{intermediate_text}\"")
    print(f"  The semantic relation of this edit is {edit1_relation} ({relation_names[edit1_relation]})")
    
    # Composition
    projected_relation_step1 = composition_table[current_relation][edit1_relation]
    print(f"  Composition: {current_relation} o {edit1_relation} = {projected_relation_step1}")
    print(f"  Projected relation after Step 1: {projected_relation_step1} ({relation_names[projected_relation_step1]})")

    # --- Step 2: Insertion ---
    edit2_relation = '^'
    print(f"\nStep 2: Edit Intermediate Text")
    print(f"  Action: Insert 'not'")
    print(f"  Final Text: \"{hypothesis}\"")
    print(f"  The semantic relation of this edit is {edit2_relation} ({relation_names[edit2_relation]})")
    
    # Composition
    final_relation = composition_table[projected_relation_step1][edit2_relation]
    print(f"  Composition: {projected_relation_step1} o {edit2_relation} = {final_relation}")
    
    # --- Final Result ---
    final_relation_name = relation_names[final_relation]
    print("-" * 25)
    print(f"The final projected relation is '{final_relation}'.")
    print(f"The name of this operator is {final_relation_name}.")

    # Output answer in specified format
    sys.stdout.flush() # Ensure all prints are out before the final answer
    print(f"<<<{final_relation_name}>>>")

solve_natural_logic_inference()