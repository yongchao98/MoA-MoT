def check_statements():
    # Assume Delfina tells the truth
    delfina_truth = True
    fletcher_truth = not delfina_truth  # Fletcher says Delfina lies
    amberly_truth = not fletcher_truth  # Amberly says Fletcher lies
    audrie_truth = amberly_truth  # Audrie says Amberly tells the truth
    lorine_truth = audrie_truth  # Lorine says Audrie tells the truth

    # Check if Lorine tells the truth
    return lorine_truth

print(check_statements())