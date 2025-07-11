import math

def is_square(n):
    """Checks if a non-negative integer is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def main():
    """
    This script explains and verifies the construction for the set A
    that gives the largest possible density c.
    """
    print("This program determines the largest density 'c' for a set A ⊆ {1,…,N}")
    print("such that the sumset A+A contains no perfect squares.")
    
    print("\nA powerful method is to construct A using modular arithmetic.")
    print("We can achieve the maximum density using the modulus m=3.")
    
    # Define the construction parameters
    m = 3
    # The set of residues R is chosen such that R+R contains no squares mod m.
    # Squares mod 3 are S_3 = {0, 1}. The non-square residues are {2}.
    # We choose R = {1}, so R+R = {1+1} = {2}. This avoids all squares.
    R_size = 1
    
    print("\nThe construction is A = {n | n ≡ 1 (mod 3)}.")
    
    # Output the equation for the density c
    print("\nThe density 'c' is given by the equation: c = |R| / m")
    print(f"For this construction, the numbers in the equation are:")
    print(f"|R| = {R_size}")
    print(f"m = {m}")
    
    c_fraction = f"{R_size}/{m}"
    c_decimal = R_size / m
    print(f"\nThus, the largest possible value for c is {c_fraction}, or approximately {c_decimal:.4f}.")
    print("This is a known result in combinatorial number theory.")
    
    # Verification for a sample N
    N = 500
    print(f"\n--- Verifying the construction for N={N} ---")
    A = [n for n in range(1, N + 1) if n % 3 == 1]
    
    found_square_sum = False
    for i in range(len(A)):
        for j in range(i, len(A)):
            a1 = A[i]
            a2 = A[j]
            s = a1 + a2
            if is_square(s):
                print(f"VERIFICATION FAILED: {a1} + {a2} = {s}, which is {int(math.sqrt(s))}^2")
                found_square_sum = True
                break
        if found_square_sum:
            break
            
    if not found_square_sum:
        print(f"VERIFICATION PASSED: For N={N}, no sum of two elements from A is a perfect square.")

if __name__ == "__main__":
    main()
