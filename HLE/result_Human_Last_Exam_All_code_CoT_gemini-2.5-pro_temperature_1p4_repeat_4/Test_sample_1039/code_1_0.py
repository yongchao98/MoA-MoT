import itertools

def calculate_inversions(p):
    """Calculates the number of inversions in a list of numbers."""
    inv_count = 0
    n = len(p)
    for i in range(n):
        for j in range(i + 1, n):
            if p[i] > p[j]:
                inv_count += 1
    return inv_count

def get_coxeter_length(w):
    """
    Calculates the Coxeter length of a signed permutation w in B_n.
    Formula: l(w) = inv(|w|) + sum_{i where w(i)<0} |w(i)|
    """
    # Part 1: Inversions of the absolute values
    abs_w = [abs(x) for x in w]
    inversions = calculate_inversions(abs_w)
    
    # Part 2: Sum of absolute values of negative entries
    neg_sum = sum(abs(x) for x in w if x < 0)
    
    return inversions + neg_sum

def main():
    """
    Calculates and prints the variance of the Coxeter length statistic
    on the hyperoctahedral group B_3.
    """
    n = 3
    
    # Step 1: Generate all elements of B_3
    # Get permutations of (1, 2, 3)
    base_perms = list(itertools.permutations(range(1, n + 1)))
    # Get all 2^3 sign combinations
    signs_list = list(itertools.product([-1, 1], repeat=n))
    
    signed_perms = []
    for p in base_perms:
        for signs in signs_list:
            # Apply signs to the permutation
            signed_perm = [p[i] * signs[i] for i in range(n)]
            signed_perms.append(signed_perm)
            
    num_elements = len(signed_perms)

    # Step 2: Calculate the length for each element
    lengths = [get_coxeter_length(w) for w in signed_perms]
    
    # Step 3: Calculate the moments (E[L] and E[L^2])
    sum_of_lengths = sum(lengths)
    mean_length = sum_of_lengths / num_elements
    
    sum_of_squared_lengths = sum(l**2 for l in lengths)
    mean_squared_length = sum_of_squared_lengths / num_elements
    
    # Step 4: Calculate the variance and display the process
    variance = mean_squared_length - mean_length**2

    print(f"The hyperoctahedral group B_3 has {num_elements} elements.")
    print("\n--- Calculating the Variance of the Coxeter Length ---")
    
    print("\n1. Mean (Expected Value) E[L]:")
    print(f"The sum of all lengths is Σl(w) = {sum_of_lengths}.")
    print(f"E[L] = (Σl(w)) / |B_3| = {sum_of_lengths} / {num_elements} = {mean_length}")

    print("\n2. Mean of Squares E[L^2]:")
    print(f"The sum of squared lengths is Σ(l(w)^2) = {sum_of_squared_lengths}.")
    print(f"E[L^2] = (Σ(l(w)^2)) / |B_3| = {sum_of_squared_lengths} / {num_elements} = {mean_squared_length:.4f}")
    
    print("\n3. Variance Var(L) = E[L^2] - (E[L])^2:")
    print(f"Var(L) = {mean_squared_length:.4f} - ({mean_length})^2")
    print(f"Var(L) = {mean_squared_length:.4f} - {mean_length**2}")
    print(f"The final variance is: {variance:.4f}")

if __name__ == "__main__":
    main()