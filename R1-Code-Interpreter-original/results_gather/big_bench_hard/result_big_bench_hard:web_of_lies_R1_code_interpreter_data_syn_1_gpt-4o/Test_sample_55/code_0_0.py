def check_truth():
    # Initial assumptions
    kandi_truth = True  # Kandi tells the truth
    alejandro_truth = kandi_truth  # Alejandro says Kandi tells the truth
    millicent_truth = not alejandro_truth  # Millicent says Alejandro lies
    shalonda_truth = not millicent_truth  # Shalonda says Millicent lies
    conception_truth = not shalonda_truth  # Conception says Shalonda lies

    return conception_truth

# Output the result
print(check_truth())