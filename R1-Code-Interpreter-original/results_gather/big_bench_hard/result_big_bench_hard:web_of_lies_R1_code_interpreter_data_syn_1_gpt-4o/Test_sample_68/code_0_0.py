def evaluate_statements():
    # Scenario 1: Andree tells the truth
    andree_truth = True
    inga_lies = True
    fletcher_lies = True
    michael_truth = True
    fidel_lies = True

    # Scenario 2: Andree lies
    andree_lies = True
    inga_truth = True
    fletcher_truth = True
    michael_lies = True
    fidel_truth = True

    # Check consistency of both scenarios
    scenario_1_consistent = inga_lies and fletcher_lies and michael_truth and fidel_lies
    scenario_2_consistent = inga_truth and fletcher_truth and michael_lies and fidel_truth

    return scenario_1_consistent, scenario_2_consistent

scenario_1, scenario_2 = evaluate_statements()
print("Scenario 1 (Andree tells the truth) is consistent:", scenario_1)
print("Scenario 2 (Andree lies) is consistent:", scenario_2)