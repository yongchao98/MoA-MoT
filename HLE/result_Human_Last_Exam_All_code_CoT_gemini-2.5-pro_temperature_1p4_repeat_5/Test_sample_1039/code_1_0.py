import itertools
from fractions import Fraction

def calculate_coxeter_length(w):
    """
    Calculates the Coxeter length of a signed permutation w.
    The formula used is: l(w) = inv(|w|) + neg(w) + cross(w)
    - inv(|w|): number of inversions in the absolute permutation.
    - neg(w): number of negative entries in w.
    - cross(w): number of pairs (i, j) with i < j such that w[i] + w[j] < 0.
    """
    n = len(w)
    abs_w = [abs(x) for x in w]

    # Term 1: inv(|w|)
    inversions = 0
    for i in range(n):
        for j in range(i + 1, n):
            if abs_w[i] > abs_w[j]:
                inversions += 1

    # Term 2: neg(w)
    neg_count = sum(1 for x in w if x < 0)

    # Term 3: cross(w)
    cross_sum_neg = 0
    for i in range(n):
        for j in range(i + 1, n):
            if w[i] + w[j] < 0:
                cross_sum_neg += 1
    
    return inversions + neg_count + cross_sum_neg

def solve_variance():
    """
    Generates all elements of the hyperoctahedral group B3, calculates their
    Coxeter lengths, and computes the variance of this distribution.
    """
    n = 3
    group_elements = []

    # Generate all base permutations of the values {1, 2, ..., n}
    base_value_perms = list(itertools.permutations(range(1, n + 1)))

    # For each permutation of values, apply all 2^n sign combinations
    for p_vals in base_value_perms:
        for signs in itertools.product([-1, 1], repeat=n):
            signed_perm = tuple(p_vals[i] * signs[i] for i in range(n))
            group_elements.append(signed_perm)

    # Calculate lengths for all elements
    lengths = [calculate_coxeter_length(w) for w in group_elements]

    # Calculate variance using fractions for exact arithmetic
    N = len(lengths)
    sum_lengths = sum(lengths)
    sum_sq_lengths = sum(l*l for l in lengths)

    mean_l = Fraction(sum_lengths, N)
    mean_sq_l = Fraction(sum_sq_lengths, N)
    variance = mean_sq_l - mean_l**2

    print("To find the variance, we use the formula: Var(X) = E[X^2] - (E[X])^2")
    print(f"\nThe group B3 has {N} elements.")
    print(f"The calculated mean of the Coxeter lengths E[X] = {sum_lengths}/{N} = {float(mean_l)}")
    print(f"The calculated mean of the squared lengths E[X^2] = {sum_sq_lengths}/{N} = {float(mean_sq_l)}")
    
    print("\nFinal variance calculation:")
    print(f"Var(X) = E[X^2] - (E[X])^2")
    print(f"Var(X) = {float(mean_sq_l)} - ({float(mean_l)})^2")
    print(f"Var(X) = {float(mean_sq_l)} - {float(mean_l**2)}")
    print(f"Var(X) = {float(variance)}")
    print(f"\nThe exact variance as a fraction is {variance.numerator}/{variance.denominator}.")

solve_variance()