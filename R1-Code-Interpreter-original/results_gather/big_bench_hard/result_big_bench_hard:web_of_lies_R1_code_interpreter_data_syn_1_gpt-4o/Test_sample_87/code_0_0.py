def check_statements():
    # Assume Teressa lies
    teressa_lies = True
    kristian_truth = teressa_lies  # Kristian says Teressa lies
    sal_truth = kristian_truth  # Sal says Kristian tells the truth
    audrie_truth = not sal_truth  # Audrie says Sal lies
    leda_truth = not audrie_truth  # Leda says Audrie lies

    return leda_truth

print(check_statements())