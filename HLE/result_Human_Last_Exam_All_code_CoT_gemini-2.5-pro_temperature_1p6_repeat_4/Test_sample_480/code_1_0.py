import sys

def solve_natural_logic_inference():
    """
    Solves the McCartney natural logic inference problem.
    This script calculates the final projected operator by simulating the
    compositional proof based on a series of edits from the premise to the hypothesis.
    """

    # MacCartney's JOIN table for composing semantic relations.
    # Rows are the outer relation (current state), columns are the inner (edit).
    JOIN_TABLE = {
        'equivalence': {'equivalence': 'equivalence', 'forward_entailment': 'forward_entailment', 'reverse_entailment': 'reverse_entailment', 'negation': 'negation', 'alternation': 'alternation', 'cover': 'cover', 'independence': 'independence'},
        'forward_entailment': {'equivalence': 'forward_entailment', 'forward_entailment': 'forward_entailment', 'reverse_entailment': 'independence', 'negation': 'independence', 'alternation': 'alternation', 'cover': 'cover', 'independence': 'independence'},
        'reverse_entailment': {'equivalence': 'reverse_entailment', 'forward_entailment': 'cover', 'reverse_entailment': 'reverse_entailment', 'negation': 'alternation', 'alternation': 'cover', 'cover': 'cover', 'independence': 'independence'},
        'negation': {'equivalence': 'negation', 'forward_entailment': 'cover', 'reverse_entailment': 'alternation', 'negation': 'equivalence', 'alternation': 'reverse_entailment', 'cover': 'forward_entailment', 'independence': 'independence'},
        'alternation': {'equivalence': 'alternation', 'forward_entailment': 'alternation', 'reverse_entailment': 'independence', 'negation': 'independence', 'alternation': 'alternation', 'cover': 'cover', 'independence': 'independence'},
        'cover': {'equivalence': 'cover', 'forward_entailment': 'forward_entailment', 'reverse_entailment': 'cover', 'negation': 'reverse_entailment', 'alternation': 'cover', 'cover': 'cover', 'independence': 'independence'},
        'independence': {'equivalence': 'independence', 'forward_entailment': 'independence', 'reverse_entailment': 'independence', 'negation': 'independence', 'alternation': 'independence', 'cover': 'independence', 'independence': 'independence'}
    }

    # MacCartney's mapping for the negation function.
    NEGATION_MAP = {
        'equivalence': 'negation',
        'forward_entailment': 'alternation',
        'reverse_entailment': 'cover',
        'negation': 'equivalence',
        'alternation': 'forward_entailment',
        'cover': 'reverse_entailment',
        'independence': 'independence'
    }

    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    print(f"Premise: \"{premise}\"")
    print(f"Hypothesis: \"{hypothesis}\"\n")
    print("Starting the compositional proof...\n")

    # Step 0: Initial state
    # The relation from a premise to itself is always equivalence.
    current_relation = 'equivalence'
    print(f"Step 0: Initial relation is {current_relation}")

    # Step 1: Handle deletion of "pop".
    # "a pop song" -> "a song". This is a generalization.
    # The lexical relation is forward_entailment (sq).
    edit1_relation = 'forward_entailment'
    previous_relation = current_relation
    current_relation = JOIN_TABLE[previous_relation][edit1_relation]
    print(f"Step 1: Deletion of 'pop' ('pop song' -> 'song'). Lexical relation is {edit1_relation}.")
    print(f"       Calculation: {previous_relation} JOIN {edit1_relation} = {current_relation}\n")


    # Step 2: Handle substitution of "Taylor Swift" with "Michael Jackson".
    # These are disjoint entities.
    # The lexical relation is alternation (alt).
    edit2_relation = 'alternation'
    previous_relation = current_relation
    current_relation = JOIN_TABLE[previous_relation][edit2_relation]
    print(f"Step 2: Substitution of 'Taylor Swift' with 'Michael Jackson'. Lexical relation is {edit2_relation}.")
    print(f"       Calculation: {previous_relation} JOIN {edit2_relation} = {current_relation}\n")

    # Step 3: Handle the insertion of "not".
    # This applies the negation function to the accumulated relation.
    previous_relation = current_relation
    current_relation = NEGATION_MAP[previous_relation]
    print(f"Step 3: Insertion of 'not'. This applies the negation function.")
    print(f"       Calculation: neg({previous_relation}) = {current_relation}\n")

    # Final Result
    final_operator_name = current_relation.replace('_', ' ').title()
    print(f"The final projected natural logic operator is: {final_operator_name}")

    # Output the final answer in the required format
    sys.stdout.write(f'<<<{final_operator_name}>>>')

solve_natural_logic_inference()