def check_statements():
    # Scenario 1: Vina tells the truth
    vina_truth = True
    helene_truth = not vina_truth
    kandi_truth = helene_truth
    jamey_truth = not kandi_truth
    ka_truth = not jamey_truth

    # Scenario 2: Vina lies
    vina_lies = not vina_truth
    helene_truth_if_vina_lies = vina_lies
    kandi_truth_if_vina_lies = helene_truth_if_vina_lies
    jamey_truth_if_vina_lies = not kandi_truth_if_vina_lies
    ka_truth_if_vina_lies = not jamey_truth_if_vina_lies

    return ka_truth, ka_truth_if_vina_lies

ka_truth, ka_truth_if_vina_lies = check_statements()
print(ka_truth, ka_truth_if_vina_lies)