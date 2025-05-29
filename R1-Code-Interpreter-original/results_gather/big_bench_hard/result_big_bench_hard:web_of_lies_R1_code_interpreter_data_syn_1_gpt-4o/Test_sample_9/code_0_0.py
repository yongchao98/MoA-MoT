def check_statements():
    # Assume Rashida tells the truth
    rashida_truth = True
    jim_truth = not rashida_truth
    alejandro_truth = not jim_truth
    tamika_truth = alejandro_truth
    fidel_truth = not tamika_truth

    # Check if Fidel tells the truth
    return fidel_truth

print(check_statements())