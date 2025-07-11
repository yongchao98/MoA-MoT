import itertools
import math
from fractions import Fraction

def solve_variance_b3():
    """
    Calculates the variance of the Coxeter length for the hyperoctahedral group B_3.
    """

    def calculate_coxeter_length(w):
        """
        Calculates the Coxeter length of a signed permutation in B_n using the formula:
        l(w) = inv(|w|) + neg(w) + npair(w)
        """
        n = len(w)
        abs_w = [abs(x) for x in w]

        # 1. inv(|w|): Number of inversions in the permutation of absolute values.
        inv_abs_w = 0
        for i in range(n):
            for j in range(i + 1, n):
                if abs_w[i] > abs_w[j]:
                    inv_abs_w += 1

        # 2. neg(w): Number of negative entries in the signed permutation.
        neg_w = sum(1 for x in w if x < 0)

        # 3. npair(w): Number of pairs (i, j) with i < j and w_i + w_j < 0.
        npair_w = 0
        for i in range(n):
            for j in range(i + 1, n):
                if w[i] + w[j] < 0:
                    npair_w += 1
                    
        return inv_abs_w + neg_w + npair_w

    # Generate all elements of B_3 and calculate their lengths
    n = 3
    base_set = list(range(1, n + 1))
    all_lengths = []

    # Iterate through permutations of {1, 2, 3}
    for p in itertools.permutations(base_set):
        # For each permutation, iterate through all 2^n sign combinations
        for signs in itertools.product([-1, 1], repeat=n):
            signed_perm = [p[i] * signs[i] for i in range(n)]
            length = calculate_coxeter_length(signed_perm)
            all_lengths.append(length)

    # Calculate the variance from the list of all lengths
    N = len(all_lengths)
    sum_of_lengths = sum(all_lengths)
    sum_of_squares = sum(x**2 for x in all_lengths)

    # Use fractions for precise calculation
    E_l = Fraction(sum_of_lengths, N)
    E_l_sq = Fraction(sum_of_squares, N)
    variance = E_l_sq - E_l**2

    # Output the explanation and the result
    print("The hyperoctahedral group B_3 has", N, "elements.")
    print("The variance is calculated as Var(l) = E[l^2] - (E[l])^2.\n")
    print(f"The sum of all lengths is: {sum_of_lengths}")
    print(f"The sum of the squares of all lengths is: {sum_of_squares}\n")

    print(f"E[l] = {sum_of_lengths}/{N} = {E_l}")
    print(f"E[l^2] = {sum_of_squares}/{N} = {E_l_sq}\n")
    
    print("Var(l) = E[l^2] - (E[l])^2")
    print(f"       = {E_l_sq} - ({E_l})^2")
    print(f"       = {E_l_sq} - {E_l**2}")
    print(f"       = {variance}")

if __name__ == '__main__':
    solve_variance_b3()