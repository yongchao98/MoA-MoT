import itertools

def coxeter_length_b_n(w):
    """
    Calculates the Coxeter length of a signed permutation w from the group B_n.
    The formula used is: l(w) = inv(|w|) + sum of absolute values of negative entries.
    """
    n = len(w)
    abs_w = [abs(x) for x in w]

    # Calculate inv(|w|), the number of inversions of the underlying permutation.
    inv_count = 0
    for i in range(n):
        for j in range(i + 1, n):
            if abs_w[i] > abs_w[j]:
                inv_count += 1

    # Calculate the sum of absolute values of all negative entries.
    neg_sum = 0
    for x in w:
        if x < 0:
            neg_sum += abs(x)

    return inv_count + neg_sum

def main():
    """
    Main function to calculate the variance of the Coxeter length for B_3.
    """
    n = 3
    
    # 1. Generate all permutations of {1, 2, ..., n}.
    base_perms = list(itertools.permutations(range(1, n + 1)))
    
    # 2. Generate all 2^n sign combinations.
    signs = list(itertools.product([-1, 1], repeat=n))
    
    all_lengths = []
    
    # 3. Generate all 2^n * n! signed permutations and calculate their lengths.
    for p in base_perms:
        for s in signs:
            w = [p[i] * s[i] for i in range(n)]
            length = coxeter_length_b_n(w)
            all_lengths.append(length)

    # 4. Calculate the components for the variance formula.
    N = len(all_lengths)
    sum_of_lengths = sum(all_lengths)
    sum_of_squares = sum(l**2 for l in all_lengths)
    
    mean = sum_of_lengths / N
    mean_of_squares = sum_of_squares / N
    variance = mean_of_squares - mean**2
    
    # 5. Output the results, showing each number in the variance equation.
    print(f"The hyperoctahedral group B_{n} has {N} elements.")
    print(f"The sum of the Coxeter lengths for all elements is: {sum_of_lengths}")
    print(f"The sum of the squares of the lengths is: {sum_of_squares}")
    print("\nThe variance is calculated using the formula: Var(X) = E[X^2] - (E[X])^2")
    
    print("\nThe final equation is:")
    print(f"Variance = ({sum_of_squares} / {N}) - ({sum_of_lengths} / {N})^2")
    print(f"         = {mean_of_squares:.4f} - ({mean:.4f})^2")
    print(f"         = {mean_of_squares:.4f} - {mean**2:.4f}")
    print(f"         = {variance:.4f}")

if __name__ == "__main__":
    main()