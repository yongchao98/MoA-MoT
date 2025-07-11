import collections

def solve_maccartney_inference():
    """
    Calculates the final projected natural logic operator for a given
    premise-hypothesis pair based on MacCartney's framework.
    """
    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    print(f"Premise:    \"{premise}\"")
    print(f"Hypothesis: \"{hypothesis}\"")
    print("-" * 30)
    print("Decomposing the inference into left-to-right atomic edits:")

    # Define the sequence of edits and their corresponding semantic relations.
    # Edits are applied from left to right as they appear in the sentences.
    # Format: (Description, Semantic Operator)
    edits = [
        ("Insert 'not' (is singing -> is not singing)", "^"), # Negation
        ("Delete 'pop' (pop song -> song)", "<"),            # Forward Entailment
        ("Substitute 'Taylor Swift' -> 'Michael Jackson'", "|") # Alternation
    ]

    # MacCartney's join table for composing semantic relations.
    operators = ['=', '<', '>', '^', '|', 'v', '#']
    join_table_data = [
        # op2 ->  =    <    >    ^    |    v    #
                ['=', '<', '>', '^', '|', 'v', '#'],  # op1 = =
                ['<', '<', '#', '>', '|', '#', '#'],  # op1 = <
                ['>', '#', '>', 'v', '#', '|', '#'],  # op1 = >
                ['^', 'v', '|', '=', '>', '<', '#'],  # op1 = ^
                ['|', '|', '#', '<', '|', '#', '#'],  # op1 = |
                ['v', '#', 'v', '>', '#', 'v', '#'],  # op1 = v
                ['#', '#', '#', '#', '#', '#', '#'],  # op1 = #
    ]
    join_table = {op1: {op2: result for op2, result in zip(operators, row)}
                  for op1, row in zip(operators, join_table_data)}

    # Map operators to their names for the final answer.
    operator_names = {
        '=': "Equivalence",
        '<': "Forward Entailment",
        '>': "Reverse Entailment",
        '^': "Negation",
        '|': "Alternation",
        'v': "Cover",
        '#': "Independence"
    }

    # Start with the Equivalence relation.
    current_relation = '='
    print(f"Initial relation is Equivalence ({current_relation})\n")

    # Sequentially compose the relations for each edit.
    for i, (description, edit_relation) in enumerate(edits):
        print(f"Step {i+1}: {description}")
        print(f"   - Semantic relation of edit: {operator_names[edit_relation]} ({edit_relation})")
        
        previous_relation = current_relation
        current_relation = join_table[previous_relation][edit_relation]
        
        print(f"   - Composition: {previous_relation} × {edit_relation} = {current_relation}")
        print(f"   - New composite relation: {operator_names[current_relation]} ({current_relation})\n")

    final_relation_name = operator_names[current_relation]
    print("-" * 30)
    print(f"The final projected natural logic operator is {current_relation}, which is named '{final_relation_name}'.")

    # Output the final answer in the specified format
    # This line is not printed to the console but is the final result for the platform.
    # print(f'<<<{final_relation_name}>>>')


solve_maccartney_inference()
<<<Independence>>>