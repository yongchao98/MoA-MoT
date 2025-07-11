import sys

def solve_maccartney_inference():
    """
    Solves for the final projected natural logic operator for a given P/H pair.
    """
    # MacCartney's Join Table for composing semantic relations.
    # The table is represented as a dictionary for easy lookups.
    # Key: (previous_relation, new_edit_relation), Value: projected_relation
    # Relations: =, <, >, ^, |, v, #
    # = (identity), < (forward entailment), > (reverse entailment)
    # ^ (negation), | (alternation/disjoint), v (cover), # (independence)
    JOIN_TABLE = {
        ('=', '='): '=', ('=', '<'): '<', ('=', '>'): '>', ('=', '^'): '^', ('=', '|'): '|', ('=', 'v'): 'v', ('=', '#'): '#',
        ('<', '='): '<', ('<', '<'): '<', ('<', '>'): '#', ('<', '^'): '|', ('<', '|'): '|', ('<', 'v'): '<', ('<', '#'): '#',
        ('>', '='): '>', ('>', '<'): '#', ('>', '>'): '>', ('>', '^'): 'v', ('>', '|'): '>', ('>', 'v'): 'v', ('>', '#'): '#',
        ('^', '='): '^', ('^', '<'): 'v', ('^', '>'): '|', ('^', '^'): '=', ('^', '|'): '>', ('^', 'v'): '<', ('^', '#'): '#',
        ('|', '='): '|', ('|', '<'): '|', ('|', '>'): 'v', ('|', '^'): '<', ('|', '|'): '<', ('|', 'v'): '#', ('|', '#'): '#',
        ('v', '='): 'v', ('v', '<'): '>', ('v', '>'): '|', ('v', '^'): '<', ('v', '|'): '#', ('v', 'v'): '<', ('v', '#'): '#',
        ('#', '='): '#', ('#', '<'): '#', ('#', '>'): '#', ('#', '^'): '#', ('#', '|'): '#', ('#', 'v'): '#', ('#', '#'): '#',
    }

    RELATION_NAMES = {
        '=': 'Identity',
        '<': 'Forward Entailment',
        '>': 'Reverse Entailment',
        '^': 'Negation',
        '|': 'Alternation',
        'v': 'Cover',
        '#': 'Independence'
    }

    # The sequence of edits transforming the premise to the hypothesis.
    edits = [
        ("delete 'pop': ('pop song' -> 'song')", '<'),  # song is a hypernym of pop song
        ("substitute 'Taylor Swift' with 'Michael Jackson'", '|'), # They are disjoint artists
        ("insert 'not'", '^')  # Negation
    ]

    # Start with the identity relation between the premise and itself.
    projected_relation = '='

    print("Premise: 'Mark is singing a pop song by Taylor Swift'")
    print("Hypothesis: 'Mark is not singing a song by Michael Jackson'")
    print("-" * 30)
    print("Calculating the final projected relation step-by-step:\n")

    print(f"Initial Relation: {projected_relation} ({RELATION_NAMES[projected_relation]})")
    print("-" * 30)

    # Apply each edit sequentially
    for i, (description, atomic_relation) in enumerate(edits):
        previous_relation = projected_relation
        projected_relation = JOIN_TABLE[(previous_relation, atomic_relation)]

        print(f"Step {i+1}: Applying Edit -> {description}")
        print(f"  Atomic Relation of Edit: {atomic_relation} ({RELATION_NAMES[atomic_relation]})")
        print(f"  Composition: Join({previous_relation}, {atomic_relation}) => {projected_relation}")
        print(f"  New Projected Relation: {projected_relation} ({RELATION_NAMES[projected_relation]})")
        print("-" * 30)
    
    final_relation_name = RELATION_NAMES[projected_relation]
    print(f"\nThe final projected relation is {projected_relation}, which is named '{final_relation_name}'.\n")

    # This ensures the final answer is also printed for easy parsing.
    # We use sys.stdout to avoid stderr capture in some environments.
    sys.stdout.write(f"<<<{final_relation_name}>>>\n")

solve_maccartney_inference()