count = 0
limit = 1000
# Sum over all possible values for x6, x5, x4
for x6 in range(limit // 6 + 1):
    for x5 in range((limit - 6*x6) // 5 + 1):
        for x4 in range((limit - 6*x6 - 5*x5) // 4 + 1):
            # Calculate x1
            x1 = limit - 4*x4 - 5*x5 - 6*x6
            # Add the number of ways to form the n_j's
            count += (x1 + 1) * (x4 + 1) * (x5 + 1)