def solve_natural_logic_inference():
    """
    Solves for the final natural logic operator based on MacCartney's framework.
    This script models the inference from Premise to Hypothesis by applying a
    series of edits and composing their semantic relations.

    Premise: "Mark is singing a pop song by Taylor Swift"
    Hypothesis: "Mark is not singing a song by Michael Jackson"
    """

    # Using MacCartney's symbols. mapping to more common symbols for clarity
    # = : equivalence
    # [ : forward entailment (subset)
    # ] : reverse entailment (superset)
    # ^ : negation (disjoint)
    # | : alternation (exhaustive, mutually exclusive)
    # v : cover (exhaustive)
    # # : independence

    # MacCartney & Manning (2009) / Angeli & Manning NatLogi join table
    JOIN_TABLE = {
        '=': {'=': '=', '[': '[', ']': ']', '^': '^', '|': '|', 'v': 'v', '#': '#'},
        '[': {'=': '[', '[': '[', ']': '#', '^': '[', '|': '[', 'v': '#', '#': '#'},
        ']': {'=': ']', '[': '#', ']': ']', '^': '^', '|': '#', 'v': 'v', '#': '#'},
        '^': {'=': '^', '[': ']', ']': '#', '^': 'v', '|': '#', 'v': 'v', '#': '#'},
        '|': {'=': '|', '[': ']', ']': '|', '^': 'v', '|': '=', 'v': 'v', '#': '#'},
        'v': {'=': 'v', '[': '#', ']': '[', '^': '#', '|': '[', 'v': 'v', '#': '#'},
        '#': {'=': '#', '[': '#', ']': '#', '^': '#', '|': '#', 'v': '#', '#': '#'},
    }

    # Projection of a relation R through a context with monotonicity M
    # ↑: upward, ↓: downward, –: non-monotone
    PROJECTION_TABLE = {
        '=': {'↑': '=', '↓': '=', '–': '#'},
        '[': {'↑': '[', '↓': ']', '–': '#'},
        ']': {'↑': ']', '↓': '[', '–': '#'},
        '^': {'↑': '^', '↓': 'v', '–': '#'},
        '|': {'↑': '|', '↓': '|', '–': '#'},
        'v': {'↑': 'v', '↓': '^', '–': '#'},
        '#': {'↑': '#', '↓': '#', '–': '#'},
    }

    RELATION_NAMES = {
        '=': "Equivalence",
        '[': "Forward Entailment",
        ']': "Reverse Entailment",
        '^': "Negation",
        '|': "Alternation",
        'v': "Cover",
        '#': "Independence"
    }

    # Start with identity relation (Premise = Premise)
    current_relation = '='
    print(f"Initial state: Premise relates to Premise with '=' ({RELATION_NAMES[current_relation]}).\n")

    # --- Edit 1: Insertion of 'not' ---
    print("Step 1: Process insertion of 'not'")
    atomic_rel_1 = '|' # A statement and its negation have an alternation relation
    print(f"  - Atomic relation is '{atomic_rel_1}' ({RELATION_NAMES[atomic_rel_1]}).")
    
    prev_relation = current_relation
    current_relation = JOIN_TABLE[prev_relation][atomic_rel_1]
    
    print(f"  - Composing relations: {prev_relation} o {atomic_rel_1} = {current_relation}")
    print(f"  - Current overall relation is now '{current_relation}' ({RELATION_NAMES[current_relation]}).\n")

    # --- Edit 2: Deletion of 'pop' ---
    print("Step 2: Process deletion of 'pop'")
    atomic_rel_2 = '[' # 'pop song' is a subset of 'song'
    monotonicity_2 = '↓' # Context is "Mark is NOT singing...", which is downward-monotone
    projected_rel_2 = PROJECTION_TABLE[atomic_rel_2][monotonicity_2]

    print(f"  - Atomic edit 'pop song' -> 'song' has relation '{atomic_rel_2}' ({RELATION_NAMES[atomic_rel_2]}).")
    print(f"  - The context 'Mark is not singing...' is downward-monotone ('{monotonicity_2}').")
    print(f"  - Projecting '{atomic_rel_2}' through '{monotonicity_2}' context yields '{projected_rel_2}'.")

    prev_relation = current_relation
    current_relation = JOIN_TABLE[prev_relation][projected_rel_2]

    print(f"  - Composing relations: {prev_relation} o {projected_rel_2} = {current_relation}")
    print(f"  - Current overall relation is now '{current_relation}' ({RELATION_NAMES[current_relation]}).\n")

    # --- Edit 3: Substitution of 'Taylor Swift' with 'Michael Jackson' ---
    print("Step 3: Process substitution 'Taylor Swift' -> 'Michael Jackson'")
    atomic_rel_3 = '^' # 'Taylor Swift' and 'Michael Jackson' are disjoint entities
    monotonicity_3 = '↓' # Context is "Mark is not singing a song by...", which is downward-monotone
    projected_rel_3 = PROJECTION_TABLE[atomic_rel_3][monotonicity_3]

    print(f"  - Atomic edit 'Taylor Swift' -> 'Michael Jackson' has relation '{atomic_rel_3}' ({RELATION_NAMES[atomic_rel_3]}).")
    print(f"  - The context 'Mark is not singing a song by...' is downward-monotone ('{monotonicity_3}').")
    print(f"  - Projecting '{atomic_rel_3}' through '{monotonicity_3}' context yields '{projected_rel_3}'.")

    prev_relation = current_relation
    final_relation = JOIN_TABLE[prev_relation][projected_rel_3]

    print(f"  - Composing relations to get the final answer:")
    # Print the final equation as requested
    print(f"    {prev_relation} o {projected_rel_3} = {final_relation}")
    
    final_relation_name = RELATION_NAMES[final_relation]
    print(f"\nThe final projected natural logic operator is '{final_relation}', which is named {final_relation_name}.")
    
    return final_relation_name

final_answer = solve_natural_logic_inference()
# The final answer is the name of the relation.
print(f"\n<<<Solution>>>\n{final_answer}\n<<<>>>")
