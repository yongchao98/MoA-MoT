def check_statements():
    # Assume Andree tells the truth
    andree_truth = True
    audrie_truth = andree_truth
    jerry_truth = not audrie_truth
    helene_truth = not jerry_truth
    bernita_truth = helene_truth

    # Check if Bernita tells the truth
    if bernita_truth:
        return "Bernita tells the truth."
    else:
        return "Bernita does not tell the truth."

print(check_statements())