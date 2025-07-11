import math

def demonstrate_square_free_sumset(N):
    """
    This function demonstrates the construction of a set A such that A+A has no squares,
    and shows its density c = |A|/N approaches 1/3.
    """
    print(f"--- Demonstration for N = {N} ---")
    
    # Step 1: Construct the set A = {n in {1, ..., N} | n = 1 (mod 3)}.
    A = {n for n in range(1, N + 1) if n % 3 == 1}
    
    # Step 2: Verify that A+A contains no perfect squares.
    # First, generate all perfect squares up to the maximum possible sum 2*N.
    max_sum = 2 * N
    squares = {k*k for k in range(1, int(math.sqrt(max_sum)) + 2)}
    
    # Check if any sum of two elements from A is a square.
    contains_square = False
    for a1 in A:
        for a2 in A:
            if (a1 + a2) in squares:
                print(f"Error: Found a square in A+A: {a1} + {a2} = {a1+a2}")
                contains_square = True
                break
        if contains_square:
            break

    if not contains_square:
        print("Verification successful: The constructed set A has the property that A+A contains no squares.")
        print("This is because any sum a1+a2 is congruent to 1+1 = 2 (mod 3), and no square is congruent to 2 (mod 3).")
    
    # Step 3: Calculate the density c = |A|/N, which relates to the equation |A| = c * N.
    size_A = len(A)
    c = size_A / N
    
    print("\nCalculating the density c:")
    print(f"The size of the set A is |A| = {size_A}")
    print(f"The size of the base set is N = {N}")
    print(f"The calculated density is c = |A| / N = {size_A} / {N} = {c:.6f}")
    
    print("\nAs N gets larger, this density approaches 1/3.")
    print(f"The value of 1/3 is approximately {1/3:.6f}")

# You can run the demonstration with a large N.
demonstrate_square_free_sumset(10000)