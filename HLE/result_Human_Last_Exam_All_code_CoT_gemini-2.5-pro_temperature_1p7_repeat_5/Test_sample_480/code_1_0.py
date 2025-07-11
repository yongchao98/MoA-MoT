import collections

# MacCartney's 7 semantic relations
# Using both names and symbols for clarity
RELATIONS = {
    "equivalence": "≡",
    "forward_entailment": "⊏", # Premise entails Hypothesis (subset)
    "reverse_entailment": "⊐", # Hypothesis entails Premise (superset)
    "negation": "^",           # Contradiction
    "alternation": "|",        # Disjoint, not exhaustive
    "cover": "v",              # Exhaustive, not disjoint
    "independence": "#"        # No relation
}
# Map names to indices for table lookup
REL_TO_IDX = {name: i for i, name in enumerate(RELATIONS.keys())}
IDX_TO_REL = {i: name for i, name in enumerate(RELATIONS.keys())}

# MacCartney's Join Table (Composition Table), represented by indices.
# Blanks in the original table are independence (#).
# JOIN[B][A] where B is the first relation and A is the second.
#      ≡  ⊏  ⊐  ^  |  v  #  (A)
JOIN = [
    [0, 1, 2, 3, 4, 5, 6], # ≡ (B)
    [1, 1, 6, 5, 6, 5, 6], # ⊏ (B)
    [2, 6, 2, 4, 4, 6, 6], # ⊐ (B)
    [3, 2, 1, 0, 5, 4, 6], # ^ (B)
    [4, 6, 2, 5, 6, 2, 6], # | (B)
    [5, 5, 6, 1, 2, 6, 6], # v (B)
    [6, 6, 6, 6, 6, 6, 6], # # (B)
]

def main():
    """
    Calculates the final projected natural logic operator for the given inference.
    """
    p = "Mark is singing a pop song by Taylor Swift"
    p_prime = "Mark is singing a song by Taylor Swift"
    h = "Mark is not singing a song by Michael Jackson"

    print("Analyzing inference:")
    print(f"  Premise (P):      '{p}'")
    print(f"  Hypothesis (H):   '{h}'")
    print("-" * 20)

    # Step 1: Relation between P and an intermediate P'
    print("Step 1: Generalize Premise P to P'")
    print(f"  P     = '{p}'")
    print(f"  P'    = '{p_prime}'")
    print("  Reasoning: 'a pop song' is a specific type of 'a song'. The predicate 'is singing' is in an upward-monotone context, so the entailment is preserved.")
    rel_p_p_prime_name = "forward_entailment"
    rel_p_p_prime_idx = REL_TO_IDX[rel_p_p_prime_name]
    print(f"  Result: Relation(P, P') is {rel_p_p_prime_name} ({RELATIONS[rel_p_p_prime_name]})")
    print("-" * 20)

    # Step 2: Relation between P' and H
    print("Step 2: Relate intermediate P' to Hypothesis H")
    print(f"  P'    = '{p_prime}'")
    print(f"  H     = '{h}'")
    print("  Reasoning: 'singing a song by TS' and 'singing a song by MJ' are mutually exclusive events (relation is Alternation, |).")
    print("  A key lemma of natural logic states that if A | B, then A ⊏ not B.")
    rel_p_prime_h_name = "forward_entailment"
    rel_p_prime_h_idx = REL_TO_IDX[rel_p_prime_h_name]
    print(f"  Result: Relation(P', H) is {rel_p_prime_h_name} ({RELATIONS[rel_p_prime_h_name]})")
    print("-" * 20)

    # Step 3: Compose the two relations
    print("Step 3: Compose the relations to find Relation(P, H)")
    final_rel_idx = JOIN[rel_p_p_prime_idx][rel_p_prime_h_idx]
    final_rel_name = IDX_TO_REL[final_rel_idx]
    final_rel_symbol = RELATIONS[final_rel_name]

    # The "numbers in the final equation" might refer to the indices of the operators being composed.
    print("  Final Composition Equation:")
    print(f"  Join(Relation(P, P'), Relation(P', H))")
    print(f"  => Join({rel_p_p_prime_name}[{rel_p_p_prime_idx}], {rel_p_prime_h_name}[{rel_p_prime_h_idx}])")
    print(f"  => {final_rel_name}[{final_rel_idx}]")
    print("-" * 20)
    
    print("Final Result:")
    print(f"The relationship from '{p}' to '{h}' is '{final_rel_name}'.")
    print(f"The final projected natural logic operator is {final_rel_name} ({final_rel_symbol}).")


if __name__ == "__main__":
    main()
