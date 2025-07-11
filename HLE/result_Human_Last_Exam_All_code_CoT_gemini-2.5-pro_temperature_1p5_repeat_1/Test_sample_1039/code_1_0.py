import itertools

def calculate_inversions(p):
    """Calculates the number of inversions in a permutation."""
    inv_count = 0
    n = len(p)
    for i in range(n):
        for j in range(i + 1, n):
            if p[i] > p[j]:
                inv_count += 1
    return inv_count

def coxeter_length_b(w):
    """
    Calculates the Coxeter length of an element w in the hyperoctahedral group B_n.
    The length is given by the formula: l(w) = inv(|w|) + sum(|w_i| for w_i < 0)
    """
    # Get the unsigned permutation |w|
    abs_w = [abs(x) for x in w]
    
    # Calculate the number of inversions in |w|
    inversions = calculate_inversions(abs_w)
    
    # Calculate the sum of absolute values of negative entries
    negative_sum = sum(abs(x) for x in w if x < 0)
    
    return inversions + negative_sum

def main():
    """
    Calculates the variance of the Coxeter length statistic on the
    hyperoctahedral group of rank 3.
    """
    n = 3
    
    # Step 1: Generate all elements of the hyperoctahedral group B_3.
    # An element is a signed permutation of {1, 2, 3}.
    base_permutations = list(itertools.permutations(range(1, n + 1)))
    sign_combinations = list(itertools.product([-1, 1], repeat=n))
    
    group_elements = []
    for p in base_permutations:
        for signs in sign_combinations:
            element = [p[i] * signs[i] for i in range(n)]
            group_elements.append(element)
            
    # Step 2: Calculate the Coxeter length for each element.
    lengths = [coxeter_length_b(w) for w in group_elements]
    
    # Step 3: Calculate the variance of the distribution of lengths.
    num_elements = len(lengths)
    
    # Calculate the mean (expected value) E[L].
    sum_of_lengths = sum(lengths)
    mean_length = sum_of_lengths / num_elements
    
    # Calculate the mean of the squares E[L^2].
    sum_of_squares = sum(l * l for l in lengths)
    mean_of_squares = sum_of_squares / num_elements
    
    # Variance Var(L) = E[L^2] - (E[L])^2
    variance = mean_of_squares - mean_length**2
    
    # Print the details of the calculation as requested.
    print(f"The rank of the hyperoctahedral group is n = {n}.")
    print(f"The size of the group B_{n} is N = {num_elements}.")
    print(f"The sum of all Coxeter lengths is ΣL = {sum_of_lengths}.")
    print(f"The mean of the Coxeter length is E[L] = ΣL / N = {sum_of_lengths}/{num_elements} = {mean_length}.")
    print(f"The sum of the squares of the lengths is ΣL² = {sum_of_squares}.")
    print(f"The mean of the squares of the lengths is E[L²] = ΣL² / N = {sum_of_squares}/{num_elements} = {mean_of_squares:.4f}.")
    print(f"The variance is Var(L) = E[L²] - (E[L])².")
    print(f"Var(L) = {mean_of_squares:.4f} - ({mean_length})² = {variance:.4f}")

if __name__ == '__main__':
    main()
