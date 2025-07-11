import itertools

def get_b_n_elements(n):
    """
    Generates all elements of the hyperoctahedral group B_n.
    An element is represented as a tuple of signed integers.
    """
    elements = []
    base_permutation = range(1, n + 1)
    
    # Iterate through all permutations of (1, 2, ..., n)
    for p in itertools.permutations(base_permutation):
        # Iterate through all 2^n sign combinations
        for signs in itertools.product([-1, 1], repeat=n):
            element = tuple(p[i] * signs[i] for i in range(n))
            elements.append(element)
    return elements

def calculate_coxeter_length(w):
    """
    Calculates the Coxeter length of a signed permutation w using the formula:
    length = inv(|w|) + sum of |w_i| for w_i < 0.
    """
    n = len(w)
    abs_w = [abs(x) for x in w]
    
    # Calculate inversions in the permutation of absolute values
    inv_count = 0
    for i in range(n):
        for j in range(i + 1, n):
            if abs_w[i] > abs_w[j]:
                inv_count += 1
                
    # Calculate sum of absolute values for negative entries
    neg_sum = 0
    for x in w:
        if x < 0:
            neg_sum += abs(x)
            
    return inv_count + neg_sum

def main():
    """
    Main function to calculate and print the variance of the Coxeter length
    for the hyperoctahedral group of rank 3.
    """
    n = 3
    
    # 1. Generate all elements of B_3
    elements = get_b_n_elements(n)
    
    # 2. Calculate the Coxeter length for each element
    lengths = [calculate_coxeter_length(w) for w in elements]
    
    # 3. Calculate the necessary sums for variance
    N = len(lengths)
    sum_of_lengths = sum(lengths)
    sum_of_squared_lengths = sum(l * l for l in lengths)
    
    # 4. Calculate mean, mean of squares, and variance
    mean_length = sum_of_lengths / N
    mean_squared_length = sum_of_squared_lengths / N
    variance = mean_squared_length - mean_length**2
    
    # 5. Print the results in a detailed equation format
    print(f"The hyperoctahedral group of rank 3, B_3, has {N} elements.")
    print(f"The sum of the Coxeter lengths for all elements is: Sum(L) = {sum_of_lengths}")
    print(f"The sum of the squares of the lengths is: Sum(L^2) = {sum_of_squared_lengths}")
    print("\nThe variance is calculated as Var(L) = E[L^2] - (E[L])^2")
    print("Step 1: Calculate the mean E[L]")
    print(f"E[L] = Sum(L) / N = {sum_of_lengths} / {N} = {mean_length}")
    print("\nStep 2: Calculate the mean of the squares E[L^2]")
    print(f"E[L^2] = Sum(L^2) / N = {sum_of_squared_lengths} / {N} = {mean_squared_length}")
    print("\nStep 3: Calculate the variance")
    print(f"Var(L) = E[L^2] - (E[L])^2 = {mean_squared_length} - ({mean_length})^2")
    print(f"Var(L) = {mean_squared_length} - {mean_length**2} = {variance}")

if __name__ == "__main__":
    main()