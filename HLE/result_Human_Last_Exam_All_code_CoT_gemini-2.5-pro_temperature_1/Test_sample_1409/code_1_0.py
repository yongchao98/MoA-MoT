import math

def solve_prism_permutations():
    """
    This script finds the number of permutations of dimensions (x, y, z)
    for a rectangular prism with volume 216, satisfying the condition
    x^(log₂(y)) = y^(log₂(z)).
    """
    
    V = 216
    
    print("The problem is to find the number of ordered triples of positive integers (x, y, z) that satisfy:")
    print(f"1. x * y * z = {V}")
    print("2. x^log₂(y) = y^log₂(z)")
    print("\nBy taking log base 2 on both sides, the second condition simplifies to:")
    print("(log₂(y)) * (log₂(x)) = (log₂(z)) * (log₂(y))")
    print("This implies that either (log₂(y) = 0) or (log₂(x) = log₂(z)).")
    print("Therefore, the condition is equivalent to (y = 1) OR (x = z).")
    
    print("\nWe use the Principle of Inclusion-Exclusion to find the total number of solutions:")
    print("Total = (Solutions where y=1) + (Solutions where x=z) - (Solutions where y=1 AND x=z)")
    
    # Case 1: y = 1
    # The equation becomes x * z = V. The number of solutions is the number of divisors of V.
    # Prime factorization of 216 = 2^3 * 3^3.
    # Number of divisors = (3+1) * (3+1) = 16.
    count_y_eq_1 = (3 + 1) * (3 + 1)
    print(f"\n1. Number of solutions where y = 1:")
    print(f"   This requires x * z = {V}. The number of solutions is the number of divisors of {V}, which is {count_y_eq_1}.")

    # Case 2: x = z
    # The equation becomes x² * y = V. We need to find the number of perfect square divisors of V.
    count_x_eq_z = 0
    solutions_x_eq_z = []
    x = 1
    while x * x <= V:
        if V % (x * x) == 0:
            y = V // (x * x)
            solutions_x_eq_z.append((x, y, x))
            count_x_eq_z += 1
        x += 1
    
    print(f"\n2. Number of solutions where x = z:")
    print(f"   This requires x² * y = {V}. The solutions (x, y, x) are:")
    for s in solutions_x_eq_z:
        print(f"   {s}")
    print(f"   There are {count_x_eq_z} such solutions.")

    # Case 3: y = 1 AND x = z (Overlap)
    # The equation becomes x² = V.
    count_overlap = 0
    sqrt_V = math.isqrt(V)
    if sqrt_V * sqrt_V == V:
        count_overlap = 1
    
    print(f"\n3. Number of solutions where y = 1 AND x = z:")
    print(f"   This requires x² = {V}. Since √{V} is not an integer, there are {count_overlap} solutions.")
    
    # Final Calculation
    total_permutations = count_y_eq_1 + count_x_eq_z - count_overlap
    
    print("\nFinal calculation:")
    print(f"Total number of permutations = {count_y_eq_1} + {count_x_eq_z} - {count_overlap}")
    print(f"Total number of permutations = {total_permutations}")

solve_prism_permutations()