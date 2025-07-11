import itertools

def solve():
    """
    Calculates the variance of the Coxeter length statistic on the
    hyperoctahedral group B_3.
    """
    n = 3

    # Generate all signed permutations for B_3. An element is a permutation of
    # {1, 2, 3} with a sign (+ or -) attached to each number.
    # The total number of elements is 2^n * n!
    base_perms = list(itertools.permutations(range(1, n + 1)))
    signs_options = list(itertools.product([-1, 1], repeat=n))
    
    signed_perms = []
    for p in base_perms:
        for s in signs_options:
            signed_perm = tuple(p[i] * s[i] for i in range(n))
            signed_perms.append(signed_perm)

    # Calculate the Coxeter length for each signed permutation.
    # The length l(w) = inv(w) + nsp(w), where inv is the number of
    # inversions in the signed sequence and nsp is the number of pairs (i, j)
    # with i <= j such that w_i + w_j < 0.
    lengths = []
    for w in signed_perms:
        # Calculate inv(w)
        inversions = 0
        for i in range(n):
            for j in range(i + 1, n):
                if w[i] > w[j]:
                    inversions += 1

        # Calculate nsp(w)
        neg_sum_pairs = 0
        for i in range(n):
            for j in range(i, n):
                if w[i] + w[j] < 0:
                    neg_sum_pairs += 1
        
        length = inversions + neg_sum_pairs
        lengths.append(length)

    # Calculate the components needed for the variance calculation.
    N = len(lengths)
    sum_of_lengths = sum(lengths)
    sum_of_squares = sum(l**2 for l in lengths)

    mean_of_lengths = sum_of_lengths / N
    mean_of_squares = sum_of_squares / N
    variance = mean_of_squares - mean_of_lengths**2
    
    # Print the explanation and the final equation with all its components.
    print("To find the variance of the Coxeter length on the hyperoctahedral group of rank 3 (B_3):")
    print("1. All 48 elements of the group were generated.")
    print("2. The Coxeter length `l` was computed for each element.")
    print("3. The variance was calculated using the formula: Var(L) = E[L^2] - (E[L])^2.")
    print("\nHere are the computed values for the formula components:")
    print(f"- Total number of elements, N = {N}")
    print(f"- Sum of all lengths, Σl = {sum_of_lengths}")
    print(f"- Sum of all squared lengths, Σl^2 = {sum_of_squares}")

    print("\nThe final variance calculation is:")
    print(f"Variance = (Σl^2 / N) - (Σl / N)^2")
    print(f"Variance = ({sum_of_squares} / {N}) - ({sum_of_lengths} / {N})^2")
    print(f"Variance = {mean_of_squares:.4f} - ({mean_of_lengths})^2")
    print(f"Variance = {mean_of_squares:.4f} - {mean_of_lengths**2:.4f}")
    print(f"Variance = {variance:.4f}")

solve()