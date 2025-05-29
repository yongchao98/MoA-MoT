def check_truth():
    # Assume Kristian tells the truth
    kristian_truth = True
    dallas_truth = kristian_truth  # Dallas says Kristian tells the truth
    osvaldo_truth = not dallas_truth  # Osvaldo says Dallas lies
    millie_truth = osvaldo_truth  # Millie says Osvaldo tells the truth
    vina_truth = millie_truth  # Vina says Millie tells the truth

    # If Kristian tells the truth, Vina should be lying
    if vina_truth:
        return "Vina tells the truth if Kristian is lying."
    else:
        return "Vina tells the truth if Kristian is lying."

print(check_truth())