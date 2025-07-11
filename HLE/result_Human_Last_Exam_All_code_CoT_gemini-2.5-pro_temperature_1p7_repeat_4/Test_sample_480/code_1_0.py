import json

def solve_entailment():
    """
    Calculates the final projected natural logic operator for the given inference.
    """
    RELATION_NAMES = {
        '=': 'Equivalence',
        '<': 'Forward Entailment',
        '>': 'Reverse Entailment',
        '^': 'Alternation',
        'v': 'Cover',
        '|': 'Independence',
        '#': 'Other (Contradiction)',
    }

    # MacCartney's Join Table for composing relations
    JOIN_TABLE = {
        '=': {'=': '=', '<': '<', '>': '>', '^': '^', 'v': 'v', '|': '|', '#': '#'},
        '<': {'=': '<', '<': '<', '>': '#', '^': '|', 'v': '<', '|': '|', '#': '#'},
        '>': {'=': '>', '<': '#', '>': '>', '^': '>', 'v': 'v', '|': 'v', '#': '#'},
        '^': {'=': '^', '<': '>', '>': '<', '^': '#', 'v': '^', '|': '>', '#': '<'},
        'v': {'=': 'v', '<': 'v', '>': '|', '^': 'v', 'v': '#', '|': '|', '#': '#'},
        '|': {'=': '|', '<': '#', '>': '<', '^': '|', 'v': '<', '|': '#', '#': '<'},
        '#': {'=': '#', '<': '#', '>': '#', '^': '#', 'v': '#', '|': '#', '#': '#'},
    }
    
    # How relations flip under negative polarity (downward monotone context)
    FLIP_RELATION = {'=': '=', '<': '>', '>': '<', '^': 'v', 'v': '^', '|': '|', '#': '#'}

    # Sequence of edits to transform Premise to Hypothesis
    # (Description, Lexical Relation, Context Polarity)
    edits = [
        ("Introduce negation 'not'", '^', '+'),
        ("Generalize 'pop song' to 'song'", '<', '-'),
        ("Substitute 'Taylor Swift' with 'Michael Jackson'", '^', '-')
    ]

    print("Premise: 'Mark is singing a pop song by Taylor Swift'")
    print("Hypothesis: 'Mark is not singing a song by Michael Jackson'\n")
    print("Starting composition chain...\n")

    current_relation = '='
    print(f"Initial Relation: = ({RELATION_NAMES[current_relation]})")

    for i, (desc, lexical_relation, polarity) in enumerate(edits):
        if polarity == '-':
            projected_relation = FLIP_RELATION[lexical_relation]
        else:
            projected_relation = lexical_relation
        
        previous_relation = current_relation
        current_relation = JOIN_TABLE[previous_relation][projected_relation]

        print(f"\nStep {i+1}: {desc}")
        print(f"  - Lexical relation is '{lexical_relation}' in context '{polarity}' -> Projects to '{projected_relation}'")
        print(f"  - Composition: join({previous_relation}, {projected_relation}) = {current_relation}")

    final_relation_name = RELATION_NAMES[current_relation]
    print("\n------------------------------------")
    print(f"Final Projected Relation: {current_relation} ({final_relation_name})")
    print("------------------------------------")

    # The final answer format as requested
    print(f"\n<<< {final_relation_name} >>>")

solve_entailment()