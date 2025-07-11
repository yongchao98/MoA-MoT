import math

def solve_for_c():
    """
    This function demonstrates the construction of a set A such that A+A
    contains no squares, to find the largest possible density c = |A|/N.
    """
    # We choose a large N to get a good approximation of the density c.
    N = 100000

    # We construct the set A using numbers congruent to 1 modulo 3.
    # For any a in A, a = 3k + 1 for some integer k.
    # The sum of any two elements a1, a2 from A is:
    # a1 + a2 = (3k+1) + (3l+1) = 3(k+l) + 2.
    # This sum is always congruent to 2 modulo 3.
    #
    # Square numbers modulo 3 can only be 0 or 1 (0*0=0, 1*1=1, 2*2=1).
    # Therefore, a sum congruent to 2 modulo 3 can never be a square.
    A = [x for x in range(1, N + 1) if x % 3 == 1]

    size_A = len(A)
    density = size_A / N

    print(f"Using N = {N}")
    print(f"We construct A = {{x in {{1..N}} | x = 1 (mod 3)}}.")
    print(f"The size of this set is |A| = {size_A}.")
    print(f"The density is |A|/N = {density:.6f}.")
    print("For large N, this density approaches 1/3.")
    print("This is known to be the largest possible density.")
    
    print("\nThe problem asks for the largest number c, which is 1/3.")
    
    # Output the final equation as requested
    numerator = 1
    denominator = 3
    value = numerator / denominator
    print(f"Final equation: {numerator} / {denominator} = {value}")

solve_for_c()