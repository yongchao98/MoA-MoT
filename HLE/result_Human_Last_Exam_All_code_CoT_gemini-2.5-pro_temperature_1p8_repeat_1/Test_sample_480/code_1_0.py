def solve_natural_logic_inference():
    """
    Solves the McCartney natural logic inference problem by chaining relations.
    """
    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"
    intermediate_hypothesis = "Mark is singing a song by Michael Jackson"

    # Define the relations
    # R1 is the relation between Premise and the Intermediate Hypothesis.
    # 'singing a pop song by TS' vs 'singing a song by MJ'. These are disjoint activities.
    relation_1_symbol = '|'
    relation_1_name = 'Alternation'

    # R2 is the relation between the Intermediate Hypothesis and the final Hypothesis.
    # One is the direct negation of the other.
    relation_2_symbol = '!'
    relation_2_name = 'Negation'

    # The chaining table C(R1, R2) gives the final relation.
    # We look up C(|, !).
    # From MacCartney's chaining table (also known as composition table):
    # C(|, !) = '<'
    final_relation_symbol = '<'
    final_relation_name = 'Forward Entailment'

    print(f"To find the relation between P and H, we use an intermediate step:")
    print(f"  P: '{premise}'")
    print(f"  H': '{intermediate_hypothesis}'")
    print(f"  H: '{hypothesis}'")
    print("-" * 20)
    print(f"Step 1: Find the relation between P and H'.")
    print(f"  The relation is '{relation_1_name}', symbolized as '{relation_1_symbol}'.")
    print("-" * 20)
    print(f"Step 2: Find the relation between H' and H.")
    print(f"  The relation is '{relation_2_name}', symbolized as '{relation_2_symbol}'.")
    print("-" * 20)
    print(f"Step 3: Chain the relations to find the relation between P and H.")
    print(f"  Final Relation = chain({relation_1_name}, {relation_2_name})")
    print(f"  Final Equation: {relation_1_symbol} o {relation_2_symbol} = {final_relation_symbol}")
    print("-" * 20)
    print(f"The final projected natural logic operator is '{final_relation_name}'.")

solve_natural_logic_inference()
<<<Forward Entailment>>>