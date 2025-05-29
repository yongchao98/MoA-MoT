def check_statements():
    # Assume Audrie lies
    audrie_truth = False
    yoland_truth = not audrie_truth
    shaunda_truth = yoland_truth
    raymond_truth = not shaunda_truth
    maybelle_truth = not raymond_truth

    # Check if the assumption leads to a consistent truth
    if maybelle_truth:
        return "Maybelle tells the truth."
    else:
        return "Maybelle does not tell the truth."

# Execute the function and print the result
print(check_statements())