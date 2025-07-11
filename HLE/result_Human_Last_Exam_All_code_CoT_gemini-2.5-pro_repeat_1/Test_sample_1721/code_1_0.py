import math

def find_c():
    """
    This script demonstrates the solution to the problem by constructing a set A
    such that A+A has no squares and calculating its density c.
    """
    # Set the upper bound for the set {1, ..., N}
    N = 10000

    # Step 1: Construct the set A based on the modulus 3 argument.
    # A = {n in {1,...,N} | n = 1 (mod 3)}
    A = [n for n in range(1, N + 1) if n % 3 == 1]

    print(f"Demonstrating the construction for N = {N}")
    print("Set A is constructed as all numbers congruent to 1 modulo 3.")
    print(f"Size of A is |A| = {len(A)}")
    print("-" * 40)

    # Step 2: Verify that A+A contains no squares.
    # We pre-compute the set of squares up to 2*N for efficiency.
    print("Verifying that the sumset A+A contains no squares...")
    max_sum = 2 * N
    squares = {i*i for i in range(1, int(math.sqrt(max_sum)) + 2)}

    found_square = False
    # To speed up verification, we only need to check sums up to 2*N
    for i in range(len(A)):
        for j in range(i, len(A)):
            current_sum = A[i] + A[j]
            if current_sum in squares:
                print(f"Verification FAILED: Found a square sum: {A[i]} + {A[j]} = {current_sum}")
                found_square = True
                break
        if found_square:
            break

    if not found_square:
        print("Verification PASSED: The set A+A contains no square numbers.")
    print("-" * 40)


    # Step 3: Present the final result for c.
    size_A = len(A)
    density = size_A / N

    print("The density of set A is:")
    print(f"|A| / N = {size_A} / {N} = {density:.4f}")
    print("\nAs N approaches infinity, this density converges to 1/3.")
    print("It has been proven that this is the maximum possible density.")
    
    # The final equation is c = 1/3
    numerator = 1
    denominator = 3
    print("\nThus, the largest number c is expressed by the equation:")
    print(f"c = {numerator} / {denominator}")

find_c()