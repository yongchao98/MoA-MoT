import sys

# Define the names for the 7 semantic relations in MacCartney's logic
RELATION_NAMES = {
    '=': 'Equivalence',
    '[': 'Forward Entailment',
    ']': 'Reverse Entailment',
    '^': 'Negation',
    '|': 'Alternation',
    '#': 'Independence',
    '~': 'Cover'
}

# MacCartney's join table for composing relations (R1 ○ R2)
# Rows are R2, Columns are R1
JOIN_TABLE = {
    # R1
    '=': {'=': '=', '[': '[', ']': ']', '^': '^', '|': '|', '#': '#', '~': '~'},
    '[': {'=': '[', '[': '[', ']': '#', '^': '|', '|': '|', '#': '#', '~': '#'},
    ']': {'=': ']', '[': '#', ']': ']', '^': '~', '|': '#', '#': '#', '~': '~'},
    '^': {'=': '^', '[': '~', ']': '|', '^': '=', '|': ']', '#': '#', '~': '['},
    '|': {'=': '|', '[': '#', ']': '[', '^': '~', '|': '=', '#': '#', '~': ']'},
    '#': {'=': '#', '[': '#', ']': '[', '^': '#', '|': '[', '#': '#', '~': '['},
    '~': {'=': '~', '[': ']', ']': '~', '^': '[', '|': '#', '#': '[', '~': '#'}
}

def solve():
    """
    Solves the entailment problem using MacCartney's natural logic framework.
    """
    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    # Sequence of edits and their corresponding relations
    edits = [
        {"original": "a pop song", "modified": "a song", "relation": '['},
        {"original": "Taylor Swift", "modified": "Michael Jackson", "relation": '|'},
        {"original": "is singing", "modified": "is not singing", "relation": '^'}
    ]

    # The initial relation is always equivalence
    current_relation = '='
    
    print(f"Premise: \"{premise}\"")
    print(f"Hypothesis: \"{hypothesis}\"")
    print("\nCalculating the final projected operator step-by-step:")
    print("-----------------------------------------------------")
    print(f"Initial Relation (Premise -> Premise): {current_relation} ({RELATION_NAMES[current_relation]})")

    # Compose the relations for each edit from left to right
    for i, edit in enumerate(edits):
        previous_relation = current_relation
        edit_relation = edit["relation"]
        current_relation = JOIN_TABLE[previous_relation][edit_relation]
        
        print(f"\nStep {i+1}: Edit \"{edit['original']}\" -> \"{edit['modified']}\"")
        print(f"  - Lexical Relation: {edit_relation} ({RELATION_NAMES[edit_relation]})")
        print(f"  - Composition: {previous_relation} ○ {edit_relation} = {current_relation}")
        print(f"  - Current Projected Relation: {current_relation} ({RELATION_NAMES[current_relation]})")
        
    print("-----------------------------------------------------")
    print(f"Final Projected Operator Symbol: {current_relation}")
    print(f"Final Projected Operator Name: {RELATION_NAMES[current_relation]}")

solve()
<<<Forward Entailment>>>