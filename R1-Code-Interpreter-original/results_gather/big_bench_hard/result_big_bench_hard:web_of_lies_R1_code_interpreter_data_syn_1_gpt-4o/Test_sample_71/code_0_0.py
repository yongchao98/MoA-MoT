def check_statements():
    # Assume Kandi tells the truth
    kandi_truth = True
    crista_truth = not kandi_truth
    ryan_truth = crista_truth
    gwenn_truth = ryan_truth
    delfina_truth = gwenn_truth

    # Check if Delfina tells the truth
    if delfina_truth:
        return "Delfina tells the truth."
    else:
        return "Delfina does not tell the truth."

print(check_statements())