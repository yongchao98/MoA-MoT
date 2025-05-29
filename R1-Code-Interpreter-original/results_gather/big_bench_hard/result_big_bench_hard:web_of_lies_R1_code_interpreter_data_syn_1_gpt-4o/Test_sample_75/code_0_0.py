def check_truth():
    # Initial assumptions
    leda_truth = True  # Leda tells the truth
    shalonda_truth = not leda_truth  # Shalonda says Leda lies
    shenna_truth = shalonda_truth  # Shenna says Shalonda tells the truth
    christie_truth = not shenna_truth  # Christie says Shenna lies
    crista_truth = christie_truth  # Crista says Christie tells the truth

    return crista_truth

print(check_truth())