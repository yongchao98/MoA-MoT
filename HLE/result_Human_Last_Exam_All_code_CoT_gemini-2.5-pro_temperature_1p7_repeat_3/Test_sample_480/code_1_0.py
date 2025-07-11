import sys

def solve():
    """
    Calculates the final projected natural logic operator for the given
    premise-hypothesis pair using MacCartney's compositional framework.
    """

    # The 7 relations and their names
    RELATIONS = {
        '=': 'equivalence',
        'sqsubset': 'forward entailment',
        'sqsupset': 'reverse entailment',
        '^': 'negation',
        '|': 'alternation',
        '#': 'cover',
        '~': 'independence',
    }

    # The join table for composing relations, JOIN(B, A)
    # B is the existing relation, A is the new relation from the edit.
    # Source: MacCartney's dissertation (2009), p. 217.
    JOIN_TABLE = {
        '=': {'=': '=', 'sqsubset': 'sqsubset', 'sqsupset': 'sqsupset', '^': '^', '|': '|', '#': '#', '~': '~'},
        'sqsubset': {'=': 'sqsubset', 'sqsubset': 'sqsubset', 'sqsupset': '~', '^': '|', '|': '|', '#': '#', '~': '~'},
        'sqsupset': {'=': 'sqsupset', 'sqsubset': '~', 'sqsupset': 'sqsupset', '^': '#', '|': '~', '#': '#', '~': '~'},
        '^': {'=': '^', 'sqsubset': 'sqsupset', 'sqsupset': 'sqsubset', '^': '=', '|': '#', '#': '|', '~': '~'},
        '|': {'=': '|', 'sqsubset': '~', 'sqsupset': '|', '^': 'sqsupset', '|': '|', '#': '=', '~': '~'},
        '#': {'=': '#', 'sqsubset': '#', 'sqsupset': '~', '^': 'sqsubset', '|': '=', '#': '#', '~': '~'},
        '~': {'=': '~', 'sqsubset': '~', 'sqsupset': '~', '^': '~', '|': '~', '#': '~', '~': '~'}
    }

    def project(relation, monotonicity):
        """Projects an atomic relation through a semantic context."""
        if monotonicity == 'upward':
            return relation
        if monotonicity == 'downward':
            if relation == 'sqsubset': return 'sqsupset'
            if relation == 'sqsupset': return 'sqsubset'
            if relation == '|': return '#'
            if relation == '#': return '|'
            return relation  # =, ^, ~ are self-dual
        return '~'  # non-monotone context as a default

    # ---- Inference Steps ----
    
    # Initial state: Premise vs. Premise
    current_relation = '='
    print("Premise: 'Mark is singing a pop song by Taylor Swift'")
    print("Hypothesis: 'Mark is not singing a song by Michael Jackson'")
    print("-" * 20)
    print("Step 0: Initialize with Premise vs. Premise.")
    print(f"Starting relation: {RELATIONS[current_relation]} ({current_relation})")
    print("-" * 20)

    # Edit 1: 'is singing' -> 'is not singing'
    atomic_rel_1 = '^'
    context_1 = 'upward'
    projected_rel_1 = project(atomic_rel_1, context_1)
    new_relation_1 = JOIN_TABLE[current_relation][projected_rel_1]

    print("Step 1: Process edit 'is singing' -> 'is not singing'.")
    print(f"  - Atomic relation is {RELATIONS[atomic_rel_1]} ({atomic_rel_1}).")
    print(f"  - Context 'Mark ___ ...' is '{context_1}' monotone.")
    print(f"  - Projected relation is {RELATIONS[projected_rel_1]} ({projected_rel_1}).")
    print(f"  - Compose: join({current_relation}, {projected_rel_1}) = {new_relation_1}")
    current_relation = new_relation_1
    print(f"Relation after Step 1: {RELATIONS[current_relation]} ({current_relation})")
    print("-" * 20)

    # Edit 2: 'a pop song by Taylor Swift' -> 'a song by Michael Jackson'
    atomic_rel_2 = '|'
    context_2 = 'downward'
    projected_rel_2 = project(atomic_rel_2, context_2)
    new_relation_2 = JOIN_TABLE[current_relation][projected_rel_2]

    print("Step 2: Process edit '... by Taylor Swift' -> '... by Michael Jackson'.")
    print(f"  - Atomic relation between phrases is {RELATIONS[atomic_rel_2]} ({atomic_rel_2}).")
    print(f"  - Context 'Mark is not singing ___' is '{context_2}' monotone.")
    print(f"  - Projected relation is {RELATIONS[projected_rel_2]} ({projected_rel_2}).")
    print(f"  - Compose: join({current_relation}, {projected_rel_2}) = {new_relation_2}")
    current_relation = new_relation_2
    print(f"Relation after Step 2: {RELATIONS[current_relation]} ({current_relation})")
    print("-" * 20)

    final_relation_name = RELATIONS[current_relation]
    print(f"The final projected natural logic operator is '{final_relation_name}'.")

solve()