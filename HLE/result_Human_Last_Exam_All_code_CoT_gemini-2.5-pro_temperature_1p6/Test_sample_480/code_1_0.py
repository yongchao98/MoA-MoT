import sys

def solve_entailment_puzzle():
    """
    Calculates the final projected natural logic operator for a given premise and hypothesis
    by simulating the compositional inference process from MacCartney's framework.
    """
    
    # Define the 7 semantic relations and their names
    relations = ["=", "⊂", "⊃", "!", "|", "#", "~"]
    RELATION_NAMES = {
        "=": "Equivalence",
        "⊂": "Forward Entailment",
        "⊃": "Reverse Entailment",
        "!": "Negation",
        "|": "Alternation",
        "#": "Independence",
        "~": "Cover",
    }
    
    # Map relation symbols to indices for table lookup
    REL_TO_IDX = {rel: i for i, rel in enumerate(relations)}

    # MacCartney's Join Table (B ∘ A) from his 2009 PhD thesis (Table 4.6)
    # This table defines the composition of two relations, R(x,z) = join(R(x,y), R(y,z)).
    # We represent it as `JOIN_TABLE[index(B)][index(A)]`.
    # Rows are B (second relation), Columns are A (first relation).
    JOIN_TABLE = [
        # A ->   =    ⊂    ⊃    !    |    #    ~
        ['=', '⊂', '⊃', '!', '|', '#', '~'],  # B = = (Equivalence)
        ['⊂', '⊂', '#', '|', '|', '#', '~'],  # B = ⊂ (Forward Entailment)
        ['⊃', '#', '⊃', '!', '#', '#', '~'],  # B = ⊃ (Reverse Entailment)
        ['!', '⊃', '⊂', '=', '⊃', '⊂', '='],  # B = ! (Negation)
        ['|', '|', '~', '⊂', '~', '⊂', '~'],  # B = | (Alternation)
        ['#', '#', '⊃', '|', '⊃', '#', '~'],  # B = # (Independence)
        ['~', '~', '⊃', '=', '⊃', '#', '='],  # B = ~ (Cover)
    ]

    def join(rel_A, rel_B):
        """Composes two relations using the JOIN_TABLE."""
        idx_A = REL_TO_IDX[rel_A]
        idx_B = REL_TO_IDX[rel_B]
        return JOIN_TABLE[idx_B][idx_A]

    # Problem setup
    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"
    
    # Sequence of edits and their corresponding effective relations
    edits = [
        {"desc": "Step 1: Negate the main verb ('is' -> 'is not').", "relation": "!", "intermediate": "Mark is not singing a pop song by Taylor Swift"},
        {"desc": "Step 2: Generalize by deleting 'pop'. (Lexical: ⊂, Context: downward) -> Effective: ⊃", "relation": "⊃", "intermediate": "Mark is not singing a song by Taylor Swift"},
        {"desc": "Step 3: Substitute disjoint items ('Taylor Swift' -> 'Michael Jackson'). (Lexical: |, Context: downward) -> Effective: |", "relation": "|", "intermediate": "Mark is not singing a song by Michael Jackson"},
    ]
    
    print("--- Calculating Final Entailment Operator ---")
    print(f"Premise (P): \"{premise}\"")
    print(f"Hypothesis (H): \"{hypothesis}\"\n")
    
    composed_relation = "="  # Start with equivalence relation P = P

    print("Composing relations step-by-step:")
    current_relation_desc = f"R(P, P) = {composed_relation} ({RELATION_NAMES[composed_relation]})"
    print(f"Initial: {current_relation_desc}\n")

    for i, edit in enumerate(edits):
        print(edit["desc"])
        
        rel_A = composed_relation
        rel_B = edit["relation"]
        composed_relation = join(rel_A, rel_B)

        print(f"  Joining: {RELATION_NAMES[rel_A]} ({rel_A}) ∘ {RELATION_NAMES[rel_B]} ({rel_B})")
        print(f"  Result of composition: {RELATION_NAMES[composed_relation]} ({composed_relation})\n")
    
    final_operator_symbol = composed_relation
    final_operator_name = RELATION_NAMES[final_operator_symbol]

    print("--- Final Result ---")
    print("The final relationship between the Premise and Hypothesis is:")
    # The user asked to output each element of the final equation
    print(f"Premise {final_operator_symbol} Hypothesis")
    
    print(f"\nThe name of the final projected natural logic operator is: {final_operator_name}")

    # Return the answer in the required format
    sys.stdout.write(f"\n<<<{final_operator_name}>>>")


solve_entailment_puzzle()