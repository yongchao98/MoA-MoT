def check_statements():
    # Assume Tamika tells the truth
    tamika_truth = True
    sima_truth = tamika_truth  # Sima says Tamika tells the truth
    rashida_truth = not sima_truth  # Rashida says Sima lies
    delfina_truth = not rashida_truth  # Delfina says Rashida lies
    jerry_truth = not delfina_truth  # Jerry says Delfina lies

    # Check if Jerry tells the truth
    return jerry_truth

print(check_statements())