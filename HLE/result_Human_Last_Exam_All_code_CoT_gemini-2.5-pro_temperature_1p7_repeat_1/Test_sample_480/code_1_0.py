def solve_mccartney_inference():
    """
    Analyzes the entailment relation for the given premise and hypothesis
    based on MacCartney's natural logic framework, printing each step.
    """
    premise = "Mark is singing a pop song by Taylor Swift"
    hypothesis = "Mark is not singing a song by Michael Jackson"
    
    print("Goal: Find the entailment relation from Premise (P) to Hypothesis (H).")
    print(f"P: \"{premise}\"")
    print(f"H: \"{hypothesis}\"")
    print("-" * 30)

    print("Step 1: Analyze the positive form of the hypothesis (H_pos).")
    print("To handle the negation cleanly, we first relate P to H_pos = ¬H.")
    h_positive = "Mark is singing a song by Michael Jackson"
    print(f"H_pos: \"{h_positive}\"")
    print("-" * 30)

    print("Step 2: Identify edits from P to H_pos and their lexical relations.")
    
    # Edit 1: From 'pop song' to 'song'
    edit1_p, edit1_h = "a pop song", "a song"
    rel1_name, rel1_sym = "Forward Entailment", "⊏"
    print(f"  - Edit 1: Substitute '{edit1_p}' with '{edit1_h}'.")
    print(f"    - Relation: '{edit1_p}' {rel1_sym} '{edit1_h}' ({rel1_name}).")
    
    # Edit 2: From 'Taylor Swift' to 'Michael Jackson'
    edit2_p, edit2_h = "Taylor Swift", "Michael Jackson"
    rel2_name, rel2_sym = "Alternation", "#"
    print(f"  - Edit 2: Substitute '{edit2_p}' with '{edit2_h}'.")
    print(f"    - Relation: '{edit2_p}' {rel2_sym} '{edit2_h}' ({rel2_name}, as they are mutually exclusive).")
    print("-" * 30)

    print("Step 3: Project and compose the relations sequentially.")
    print("Since the context in P ('Mark is singing...') is positive (upward monotone), the projected relations are the same as the lexical ones.")
    
    # Composition
    initial_rel_sym = "≡"
    # The equation for the first composition
    print("\nEquation 1: Compose initial relation with the first projected relation.")
    print(f"  {initial_rel_sym} ∘ {rel1_sym} = {rel1_sym}")
    intermediate_relation_sym = rel1_sym # The result is ⊏
    
    # The equation for the second composition
    final_phpos_rel_sym = "#"
    print("\nEquation 2: Compose the result with the second projected relation.")
    print(f"  Derivation: A ⊏ B (A ⇒ B) and B # C (B ⇒ ¬C), implies A ⇒ ¬C. This is A # C.")
    print(f"  {intermediate_relation_sym} ∘ {rel2_sym} = {final_phpos_rel_sym}")
    print("-" * 30)
    
    print("Step 4: State the relation between P and H_pos.")
    print(f"The final composed relation between P and H_pos is Alternation ({final_phpos_rel_sym}).")
    print("This means: P entails the negation of H_pos.")
    print(f"Final Equation for P vs H_pos: P ⇒ ¬(H_pos)")
    print("-" * 30)
    
    print("Step 5: Convert the result to the final relation between P and H.")
    print("We know that H = ¬(H_pos). Substituting this into our equation:")
    print("Final Equation for P vs H: P ⇒ ¬(¬H)")
    print("This simplifies to the final relationship: P ⇒ H")
    
    final_operator_name = "Forward Entailment"
    final_operator_symbol = "⊏"
    print("\nThe relationship 'P entails H' (P ⇒ H) corresponds to the operator:")
    print(f"--> {final_operator_name} ({final_operator_symbol})")


if __name__ == "__main__":
    solve_mccartney_inference()