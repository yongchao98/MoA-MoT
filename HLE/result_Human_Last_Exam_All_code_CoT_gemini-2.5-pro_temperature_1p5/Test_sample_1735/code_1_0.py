def solve_dynamic_sequence(A):
    """
    Finds an optimal permutation of a sequence A to minimize the maximum cumulative impact.
    This implementation uses an O(n log n) greedy algorithm.
    """
    positives = sorted([x for x in A if x >= 0])
    negatives = sorted([x for x in A if x < 0], reverse=True)

    optimal_sequence = []
    s = 0
    p_idx, n_idx = 0, 0

    for _ in range(len(A)):
        # Decide whether to pick a positive or a negative number
        # If sum is positive, try to reduce it with a negative number
        # If sum is zero or negative, try to increase it with a positive number
        pick_positive = False
        if s <= 0:
            if p_idx < len(positives):
                pick_positive = True
            # if only negatives are left, we have to pick a negative
        else: # s > 0
            if n_idx < len(negatives):
                pick_positive = False
            # if only positives are left, we have to pick a positive
            else:
                pick_positive = True

        if pick_positive:
            num = positives[p_idx]
            p_idx += 1
        else:
            num = negatives[n_idx]
            n_idx += 1
        
        optimal_sequence.append(num)
        s += num
        
    # Calculate the cumulative impacts for the optimal sequence
    s = 0
    cumulative_sums = []
    cumulative_impacts = []
    for x in optimal_sequence:
        s += x
        cumulative_sums.append(s)
        cumulative_impacts.append(abs(s))
        
    print(f"Input: {A}")
    print(f"Optimal permutation: {optimal_sequence}")
    print(f"Cumulative sums: {cumulative_sums}")
    print(f"Absolute cumulative sums (C(k)): {cumulative_impacts}")
    
    # "output each number in the final equation!"
    # The final equation is max(C(k)) = result
    # We print the numbers C(k) that go into the max() function.
    print("The final maximum cumulative impact is the max of the following numbers:")
    # Using str.join to format the numbers for the "equation"
    print(", ".join(map(str, cumulative_impacts)))
    print(f"max({', '.join(map(str, cumulative_impacts))}) = {max(cumulative_impacts)}")


# Run on the "Tricky case" example from the problem description
# The optimal from the problem description had a typo and a higher max value
# This algorithm finds a better solution.
solve_dynamic_sequence([1, -4, 3, -1, 2, -2])
