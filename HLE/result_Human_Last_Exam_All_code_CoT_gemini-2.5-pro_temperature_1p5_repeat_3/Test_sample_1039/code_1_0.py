import itertools

def calculate_variance_of_coxeter_length():
    """
    Generates all elements of the hyperoctahedral group B_3, calculates
    the Coxeter length for each, and then computes the variance of this statistic.
    """
    # Step 1: Generate all elements of B_3
    elements = []
    base_permutations = list(itertools.permutations([1, 2, 3]))
    # Generate all 2^3 = 8 sign combinations
    sign_combinations = list(itertools.product([-1, 1], repeat=3))

    for p in base_permutations:
        for s in sign_combinations:
            element = tuple(p[i] * s[i] for i in range(3))
            elements.append(element)
    
    # Step 2 & 3: Calculate lengths and their sums
    lengths = []
    sum_of_lengths = 0
    sum_of_squared_lengths = 0
    n = len(elements)

    for w in elements:
        # Calculate inv(w): number of pairs (i, j) with i < j and w[i] > w[j]
        inversions = 0
        for i in range(len(w)):
            for j in range(i + 1, len(w)):
                if w[i] > w[j]:
                    inversions += 1
        
        # Calculate nsum(w): sum of absolute values of negative entries
        neg_sum = sum(abs(x) for x in w if x < 0)
        
        length = inversions + neg_sum
        lengths.append(length)
        sum_of_lengths += length
        sum_of_squared_lengths += length**2

    # Step 4 & 5: Compute mean, mean of squares, and variance
    if n == 0:
        mean_length = 0
        mean_squared_length = 0
    else:
        mean_length = sum_of_lengths / n
        mean_squared_length = sum_of_squared_lengths / n
    
    variance = mean_squared_length - mean_length**2

    # Print the results as a step-by-step equation
    print(f"The hyperoctahedral group B_3 has {n} elements.")
    print(f"The sum of all Coxeter lengths is: {sum_of_lengths}")
    print(f"The sum of all squared Coxeter lengths is: {sum_of_squared_lengths}")
    print("\nCalculating the Variance:")
    print("Variance = E[length^2] - (E[length])^2\n")
    print(f"E[length] = {sum_of_lengths} / {n} = {mean_length}")
    print(f"E[length^2] = {sum_of_squared_lengths} / {n} = {mean_squared_length}")
    print(f"\nVariance = {mean_squared_length} - ({mean_length})^2")
    print(f"Variance = {mean_squared_length} - {mean_length**2}")
    print(f"Variance = {variance}")

if __name__ == '__main__':
    calculate_variance_of_coxeter_length()