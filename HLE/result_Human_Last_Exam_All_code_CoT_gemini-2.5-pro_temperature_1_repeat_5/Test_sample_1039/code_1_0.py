import itertools

def calculate_b3_variance():
    """
    Calculates the variance of the Coxeter length for the hyperoctahedral group B_3.
    """
    n = 3
    base_permutations = list(itertools.permutations(range(1, n + 1)))
    signed_permutations = []

    # 1. Generate all 48 signed permutations of B_3
    for p in base_permutations:
        # There are 2^n sign combinations for each permutation of absolute values
        for signs in itertools.product([-1, 1], repeat=n):
            signed_perm = [p[i] * signs[i] for i in range(n)]
            signed_permutations.append(signed_perm)

    lengths = []
    # 2. Calculate the Coxeter length for each element
    for w in signed_permutations:
        # l(w) = inv(w) + nsum(w)
        
        # Calculate inv(w): number of pairs (i, j) with i < j and w[i] > w[j]
        inversions = 0
        for i in range(n):
            for j in range(i + 1, n):
                if w[i] > w[j]:
                    inversions += 1
        
        # Calculate nsum(w): sum of absolute values of negative entries
        nsum = 0
        for x in w:
            if x < 0:
                nsum += abs(x)
        
        length = inversions + nsum
        lengths.append(length)

    # 3. Calculate the variance of the list of lengths
    count = len(lengths)
    
    # E[L] = mean of the lengths
    sum_of_lengths = sum(lengths)
    mean = sum_of_lengths / count
    
    # E[L^2] = mean of the squares of the lengths
    sum_of_squares = sum([x**2 for x in lengths])
    mean_sq = sum_of_squares / count
    
    # Var(L) = E[L^2] - (E[L])^2
    variance = mean_sq - mean**2

    print("The variance is calculated using the formula: Var(L) = E[L^2] - (E[L])^2")
    print(f"The mean of the square of the Coxeter length is E[L^2] = {sum_of_squares} / {count} = {mean_sq}")
    print(f"The mean of the Coxeter length is E[L] = {sum_of_lengths} / {count} = {mean}")
    print(f"The final variance is {mean_sq} - ({mean})^2 = {variance}")

if __name__ == '__main__':
    calculate_b3_variance()
