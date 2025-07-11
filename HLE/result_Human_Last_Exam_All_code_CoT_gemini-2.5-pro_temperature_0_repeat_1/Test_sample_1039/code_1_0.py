import itertools

def calculate_b3_variance():
    """
    Calculates the variance of the Coxeter length for the hyperoctahedral group B3.
    """
    n = 3
    
    # Step 1: Generate all 48 elements of the group B3.
    # B3 is the group of signed permutations of {1, 2, 3}.
    base_permutations = list(itertools.permutations(range(1, n + 1)))
    sign_combinations = list(itertools.product([-1, 1], repeat=n))
    
    all_elements = []
    for p in base_permutations:
        for s in sign_combinations:
            element = [p[i] * s[i] for i in range(n)]
            all_elements.append(element)
            
    # Step 2: Calculate the Coxeter length for each element.
    lengths = []
    for w in all_elements:
        # The length l(w) is the sum of three components:
        # neg(w): number of negative entries
        # inv(w): number of standard inversions
        # nsp(w): number of negative sum pairs
        
        neg_w = sum(1 for x in w if x < 0)
        
        inv_w = 0
        for i in range(n):
            for j in range(i + 1, n):
                if w[i] > w[j]:
                    inv_w += 1
                    
        nsp_w = 0
        for i in range(n):
            for j in range(i + 1, n):
                if w[i] + w[j] < 0:
                    nsp_w += 1
                    
        length = neg_w + inv_w + nsp_w
        lengths.append(length)
        
    # Step 3: Compute the variance of the lengths.
    # The variance is E[L^2] - (E[L])^2.
    
    N = len(lengths)
    sum_of_lengths = sum(lengths)
    sum_of_squared_lengths = sum(l**2 for l in lengths)
    
    mean_length = sum_of_lengths / N
    variance = (sum_of_squared_lengths / N) - (mean_length ** 2)
    
    # Output the numbers used in the final variance equation.
    print("The variance is calculated using the formula: Var(L) = S2/N - (S1/N)^2")
    print(f"Number of elements (N): {N}")
    print(f"Sum of lengths (S1): {sum_of_lengths}")
    print(f"Sum of squared lengths (S2): {sum_of_squared_lengths}")
    print("\nFinal Equation:")
    print(f"Variance = {sum_of_squared_lengths}/{N} - ({sum_of_lengths}/{N})^2")
    print(f"Variance = {sum_of_squared_lengths/N} - ({mean_length})^2")
    print(f"Variance = {variance}")

if __name__ == '__main__':
    calculate_b3_variance()