import sys

def solve_entailment():
    """
    Calculates the final projected natural logic operator for the given inference
    using a simplified model of MacCartney's compositional semantics.
    """

    # Define the 7 logic relations and their names
    RELATION_NAMES = {
        '=': "Equivalence",
        '<': "Forward Entailment",
        '>': "Reverse Entailment",
        '^': "Exclusion",
        '|': "Cover",
        '#': "Independence"
    }

    # MacCartney's Join Table J(r1, r2) from his 2009 thesis
    # Rows are r1, columns are r2
    JOIN_TABLE = {
        # r2 -->
        # =    <    >    ^    |    #      (r1_rows)
        '=': {'=': '=', '<': '<', '>': '>', '^': '^', '|': '|', '#': '#'}, # =
        '<': {'=': '<', '<': '<', '>': '#', '^': '^', '|': '|', '#': '#'}, # <
        '>': {'=': '>', '<': '#', '>': '>', '^': '#', '|': '>', '#': '#'}, # >
        '^': {'=': '^', '<': '^', '>': '#', '^': '=', '|': '<', '#': '#'}, # ^
        '|': {'=': '|', '<': '#', '>': '>', '^': '<', '|': '=', '#': '>'}, # |
        '#': {'=': '#', '<': '#', '>': '#', '^': '#', '|': '#', '#': '#'}, # #
    }

    # Negation map: neg(R(x, y)) gives R(x, not y)
    # The relation between Taylor Swift and Michael Jackson is Exclusion (^).
    # 'x is TS' entails 'x is not MJ'. Thus, neg(^) should be <.
    NEGATION_MAP = {
        '=': '^',
        '<': '|',
        '>': '#',
        '^': '<', # Crucial assumption for exclusion
        '|': '<',
        '#': '>'
    }

    print("Step 1: Define the premise (P) and the positive version of the hypothesis (H').")
    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis_prime = "Mark is singing a song by Michael Jackson"
    print(f"P: \"{premise}\"")
    print(f"H': \"{hypothesis_prime}\"\n")
    
    print("Step 2: Determine lexical relations for edits transforming P to H'.")
    # Edits are processed in a stable, upward-monotone context.
    edit1_lexical = '<' # pop song -> song
    edit2_lexical = '^' # Taylor Swift -> Michael Jackson
    print(f"Edit 1 ('pop song' -> 'song'): {edit1_lexical} ({RELATION_NAMES[edit1_lexical]})")
    print(f"Edit 2 ('Taylor Swift' -> 'Michael Jackson'): {edit2_lexical} ({RELATION_NAMES[edit2_lexical]})\n")

    print("Step 3: Compose relations to find the relation between P and H'.")
    # The context "Mark is singing..." is upward monotone, so projection is identity.
    # We start with an identity relation.
    current_relation = '='
    print(f"Initial Relation R_0: {current_relation}")

    # Join relation from Edit 1
    rel_after_edit1 = JOIN_TABLE[current_relation][edit1_lexical]
    print(f"Equation for R_1: R_1 = J(R_0, R_edit1) = J({current_relation}, {edit1_lexical}) = {rel_after_edit1}")
    
    # Join relation from Edit 2
    rel_after_edit2 = JOIN_TABLE[rel_after_edit1][edit2_lexical]
    print(f"Equation for R_2: R_2 = J(R_1, R_edit2) = J({rel_after_edit1}, {edit2_lexical}) = {rel_after_edit2}")
    relation_P_H_prime = rel_after_edit2
    print(f"\nThe relation R(P, H') is: {relation_P_H_prime} ({RELATION_NAMES[relation_P_H_prime]})\n")

    print("Step 4: Apply negation to find the final relation between P and H.")
    print("The final hypothesis H is 'not H''.")
    final_relation = NEGATION_MAP[relation_P_H_prime]
    print(f"Final Equation: R(P, H) = neg(R(P, H')) = neg({relation_P_H_prime}) = {final_relation}\n")
    
    print("---")
    print("Final Answer")
    print("The final projected natural logic operator symbol is:", final_relation)
    print("The name of the operator is:", RELATION_NAMES[final_relation])

solve_entailment()