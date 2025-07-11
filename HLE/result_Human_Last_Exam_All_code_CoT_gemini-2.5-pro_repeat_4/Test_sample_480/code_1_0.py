import sys

def solve_inference():
    """
    Calculates the final projected natural logic operator for the given premise and hypothesis
    based on MacCartney's framework.
    """
    # Define relation symbols, their mapping to indices for table lookup, and their names.
    REL_TO_IDX = {
        "≡": 0, "⊂": 1, "⊃": 2, "^": 3, "|": 4, "⌣": 5, "#": 6
    }
    IDX_TO_SYM = {v: k for k, v in REL_TO_IDX.items()}
    IDX_TO_NAME = {
        0: "equivalence", 1: "forward entailment", 2: "reverse entailment",
        3: "negation", 4: "alternation", 5: "cover", 6: "independence"
    }

    # Define MacCartney's join table for projecting a current relation (row)
    # over an edit relation (column).
    # Order: ≡, ⊂, ⊃, ^, |, ⌣, #
    JOIN_TABLE = [
        # ≡,  ⊂,  ⊃,  ^,  |,  ⌣,  #  <- Edit Relation
        [ 0,  1,  2,  3,  4,  5,  6], # ≡ (Current Relation)
        [ 1,  1,  6,  4,  4,  6,  6], # ⊂
        [ 2,  6,  2,  5,  6,  5,  6], # ⊃
        [ 3,  5,  6,  0,  2,  6,  6], # ^
        [ 4,  6,  1,  5,  6,  1,  6], # |
        [ 5,  5,  6,  1,  2,  6,  6], # ⌣
        [ 6,  6,  6,  6,  6,  6,  6]  # #
    ]

    def project(current_rel_idx, edit_rel_idx):
        """Projects the current relation over the edit relation using the join table."""
        return JOIN_TABLE[current_rel_idx][edit_rel_idx]

    # Define premise, hypothesis, and the sequence of edits with their lexical relations.
    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"

    edits = [
        {
            "description": "Delete 'pop' (pop song → song)",
            "relation_symbol": "⊂" # forward entailment
        },
        {
            "description": "Substitute 'Taylor Swift' with 'Michael Jackson'",
            "relation_symbol": "|" # alternation
        },
        {
            "description": "Negate 'is singing' (is singing → is not singing)",
            "relation_symbol": "^" # negation
        }
    ]

    # Start with equivalence relation and apply edits sequentially.
    current_relation_symbol = "≡"
    print(f"Premise: \"{premise}\"")
    print(f"Hypothesis: \"{hypothesis}\"\n")
    print(f"Starting with initial relation: {current_relation_symbol} ({IDX_TO_NAME[REL_TO_IDX[current_relation_symbol]]})")
    print("-" * 40)

    for i, edit in enumerate(edits):
        print(f"Step {i+1}: {edit['description']}")
        edit_relation_symbol = edit['relation_symbol']
        print(f"Lexical relation for this edit: {edit_relation_symbol} ({IDX_TO_NAME[REL_TO_IDX[edit_relation_symbol]]})")

        current_rel_idx = REL_TO_IDX[current_relation_symbol]
        edit_rel_idx = REL_TO_IDX[edit_relation_symbol]

        new_rel_idx = project(current_rel_idx, edit_rel_idx)
        new_relation_symbol = IDX_TO_SYM[new_rel_idx]

        # Output the projection equation for this step
        print(f"Projecting: {current_relation_symbol} ⨝ {edit_relation_symbol} = {new_relation_symbol}")

        current_relation_symbol = new_relation_symbol
        print(f"New relation is: {current_relation_symbol} ({IDX_TO_NAME[REL_TO_IDX[current_relation_symbol]]})")
        print("-" * 40)

    final_relation_name = IDX_TO_NAME[REL_TO_IDX[current_relation_symbol]]
    print(f"The final projected natural logic operator is '{final_relation_name}'.")

    # Output final answer in the specified format.
    # We use sys.stdout.write to avoid adding a newline after the answer tag.
    sys.stdout.write(f"<<<{final_relation_name}>>>")

solve_inference()