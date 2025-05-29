def check_statements():
    # Assume Audrie tells the truth
    audrie_truth = True
    yoland_truth = not audrie_truth
    shaunda_truth = yoland_truth
    raymond_truth = not shaunda_truth
    maybelle_truth = not raymond_truth

    # Check if the assumption leads to a contradiction
    if maybelle_truth:
        return "Maybelle tells the truth."
    else:
        return "Maybelle does not tell the truth."

# Execute the function and print the result
print(check_statements())