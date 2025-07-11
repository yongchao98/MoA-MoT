import itertools
import math

def calculate_coxeter_length(w):
    """
    Calculates the Coxeter length of a signed permutation w in B_n.
    The formula is l(w) = inv(w) + nsp(w).
    
    Args:
        w (tuple): A signed permutation, e.g., (-1, 3, -2).
        
    Returns:
        int: The Coxeter length of w.
    """
    n = len(w)
    
    # Calculate inv(w), the number of inversions on the signed sequence.
    # An inversion is a pair (i, j) with i < j such that w[i] > w[j].
    inversions = 0
    for i in range(n):
        for j in range(i + 1, n):
            if w[i] > w[j]:
                inversions += 1
                
    # Calculate nsp(w), the number of negative sum pairs.
    # A negative sum pair is a pair (i, j) with i <= j such that w[i] + w[j] < 0.
    neg_sum_pairs = 0
    for i in range(n):
        for j in range(i, n):
            if w[i] + w[j] < 0:
                neg_sum_pairs += 1
                
    return inversions + neg_sum_pairs

def solve():
    """
    Calculates the variance of the Coxeter length for the hyperoctahedral group B_3.
    """
    n = 3
    
    # Step 1: Generate all elements of B_3.
    # Start with permutations of (1, 2, 3).
    base_permutations = list(itertools.permutations(range(1, n + 1)))
    
    # For each permutation, apply all possible sign combinations.
    signed_permutations = []
    for p in base_permutations:
        # There are 2^n possible sign combinations for a permutation of length n.
        for i in range(2**n):
            signed_p = list(p)
            # Use the bits of 'i' to determine the signs.
            for j in range(n):
                if (i >> j) & 1:
                    signed_p[j] *= -1
            signed_permutations.append(tuple(signed_p))
            
    # Step 2 & 3: Calculate the Coxeter length for each element.
    lengths = [calculate_coxeter_length(w) for w in signed_permutations]
    
    num_elements = len(lengths)
    
    # Step 4: Compute E[L] and E[L^2].
    sum_of_lengths = sum(lengths)
    mean_length = sum_of_lengths / num_elements
    
    sum_of_squared_lengths = sum(l**2 for l in lengths)
    mean_squared_length = sum_of_squared_lengths / num_elements
    
    # Step 5: Calculate the variance.
    variance = mean_squared_length - mean_length**2
    
    # Print the results of the calculation.
    print(f"The number of elements in the hyperoctahedral group B_{n} is {num_elements}.")
    print(f"The sum of the Coxeter lengths is {sum_of_lengths}.")
    print(f"The mean of the Coxeter length is E[L] = {sum_of_lengths} / {num_elements} = {mean_length}.")
    print(f"The sum of the squared Coxeter lengths is {sum_of_squared_lengths}.")
    print(f"The mean of the squared Coxeter length is E[L^2] = {sum_of_squared_lengths} / {num_elements} = {mean_squared_length}.")
    print("\nThe variance is calculated as Var(L) = E[L^2] - (E[L])^2")
    print(f"Var(L) = {mean_squared_length} - ({mean_length})^2 = {mean_squared_length} - {mean_length**2} = {variance}")

solve()
<<<8.5>>>