import sys

def solve_entailment_projection():
    """
    Solves the MacCartney natural logic projection problem.

    This function simulates the compositional inference process to find the
    final projected entailment relation between the given premise and hypothesis.
    It uses MacCartney's 7-relation model and the corresponding projection table.
    """

    # MacCartney's 7 semantic relations
    RELATIONS = {
        '=': {'name': 'equivalence', 'idx': 0},
        '<': {'name': 'forward entailment', 'idx': 1},
        '>': {'name': 'reverse entailment', 'idx': 2},
        '^': {'name': 'negation', 'idx': 3},
        '|': {'name': 'alternation', 'idx': 4},
        'v': {'name': 'cover', 'idx': 5},
        '#': {'name': 'independence', 'idx': 6},
    }
    
    # Create a reverse mapping from index to symbol and name
    IDX_TO_REL = {v['idx']: (k, v['name']) for k, v in RELATIONS.items()}

    # Dual relations for flipping in downward-monotone contexts
    DUAL_RELATIONS = {
        RELATIONS['=']['idx']: RELATIONS['=']['idx'],
        RELATIONS['<']['idx']: RELATIONS['>']['idx'],
        RELATIONS['>']['idx']: RELATIONS['<']['idx'],
        RELATIONS['^']['idx']: RELATIONS['^']['idx'],
        RELATIONS['|']['idx']: RELATIONS['v']['idx'],
        RELATIONS['v']['idx']: RELATIONS['|']['idx'],
        RELATIONS['#']['idx']: RELATIONS['#']['idx'],
    }

    # MacCartney's projection join table (from his 2009 dissertation, Figure 5.4)
    # Rows are the current projected relation (p), columns are the edit's relation (r).
    # JOIN_TABLE[p][r]
    # Indices: 0: =, 1: <, 2: >, 3: ^, 4: |, 5: v, 6: #
    JOIN_TABLE = [
        # r =  =  <  >  ^  |  v  #
        [0, 1, 2, 3, 4, 5, 6],  # p = =
        [1, 1, 6, 4, 4, 6, 6],  # p = <
        [2, 6, 2, 5, 6, 5, 6],  # p = >
        [3, 5, 4, 0, 2, 1, 6],  # p = ^
        [4, 6, 2, 1, 6, 2, 6],  # p = |
        [5, 1, 6, 2, 1, 6, 6],  # p = v
        [6, 6, 6, 6, 6, 6, 6],  # p = #
    ]

    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    print(f"Premise: \"{premise}\"")
    print(f"Hypothesis: \"{hypothesis}\"\n")

    # Start with the identity relation: Premise entails itself.
    current_relation_idx = RELATIONS['=']['idx']
    print(f"Step 0: Initial relation is {IDX_TO_REL[current_relation_idx][1]} ({IDX_TO_REL[current_relation_idx][0]}).\n")

    # --- Edit 1: 'is singing' -> 'is not singing' ---
    step = 1
    edit_desc = "'is singing' -> 'is not singing'"
    atomic_rel_idx = RELATIONS['^']['idx']
    # Context is upward-monotone, so the relation does not flip.
    effective_rel_idx = atomic_rel_idx
    prev_relation_idx = current_relation_idx
    current_relation_idx = JOIN_TABLE[prev_relation_idx][effective_rel_idx]

    print(f"Step {step}: Process edit \"{edit_desc}\"")
    print(f"  - Atomic relation is {IDX_TO_REL[atomic_rel_idx][1]} ({IDX_TO_REL[atomic_rel_idx][0]}).")
    print(f"  - Context is upward-monotone. Effective relation is unchanged.")
    print(f"  - Projecting: project({IDX_TO_REL[prev_relation_idx][1]}, {IDX_TO_REL[effective_rel_idx][1]}) = {IDX_TO_REL[current_relation_idx][1]}")
    print(f"  - Current overall relation is now: {IDX_TO_REL[current_relation_idx][1]} ({IDX_TO_REL[current_relation_idx][0]})\n")

    # --- Edit 2: 'pop song' -> 'song' ---
    step = 2
    edit_desc = "'pop song' -> 'song'"
    atomic_rel_idx = RELATIONS['<']['idx']
    # Context '...is not...' is downward-monotone, so the relation flips to its dual.
    effective_rel_idx = DUAL_RELATIONS[atomic_rel_idx]
    prev_relation_idx = current_relation_idx
    current_relation_idx = JOIN_TABLE[prev_relation_idx][effective_rel_idx]

    print(f"Step {step}: Process edit \"{edit_desc}\"")
    print(f"  - Atomic relation is {IDX_TO_REL[atomic_rel_idx][1]} ({IDX_TO_REL[atomic_rel_idx][0]}).")
    print(f"  - Context 'Mark is not...' is downward-monotone. Relation flips to its dual: {IDX_TO_REL[effective_rel_idx][1]} ({IDX_TO_REL[effective_rel_idx][0]}).")
    print(f"  - Projecting: project({IDX_TO_REL[prev_relation_idx][1]}, {IDX_TO_REL[effective_rel_idx][1]}) = {IDX_TO_REL[current_relation_idx][1]}")
    print(f"  - Current overall relation is now: {IDX_TO_REL[current_relation_idx][1]} ({IDX_TO_REL[current_relation_idx][0]})\n")

    # --- Edit 3: 'Taylor Swift' -> 'Michael Jackson' ---
    step = 3
    edit_desc = "'Taylor Swift' -> 'Michael Jackson'"
    atomic_rel_idx = RELATIONS['|']['idx']
    # Context '...is not...' is downward-monotone, so the relation flips to its dual.
    effective_rel_idx = DUAL_RELATIONS[atomic_rel_idx]
    prev_relation_idx = current_relation_idx
    current_relation_idx = JOIN_TABLE[prev_relation_idx][effective_rel_idx]
    
    final_op_name = IDX_TO_REL[current_relation_idx][1]

    print(f"Step {step}: Process edit \"{edit_desc}\"")
    print(f"  - Atomic relation is {IDX_TO_REL[atomic_rel_idx][1]} ({IDX_TO_REL[atomic_rel_idx][0]}).")
    print(f"  - Context 'Mark is not...' is downward-monotone. Relation flips to its dual: {IDX_TO_REL[effective_rel_idx][1]} ({IDX_TO_REL[effective_rel_idx][0]}).")
    print(f"  - Final Equation: project({IDX_TO_REL[prev_relation_idx][1]}, {IDX_TO_REL[effective_rel_idx][1]}) = {final_op_name}")
    print(f"  - The final projected relation is: {final_op_name} ({IDX_TO_REL[current_relation_idx][0]})\n")

    # The final answer in the specified format
    # Redirecting to stdout to ensure it's the last line for parsing
    sys.stdout.write(f"<<<{final_op_name}>>>\n")

solve_entailment_projection()