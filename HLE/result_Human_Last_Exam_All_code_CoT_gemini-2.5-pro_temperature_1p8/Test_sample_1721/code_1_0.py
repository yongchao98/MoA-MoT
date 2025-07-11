import math

def solve_square_sum_problem(N):
    """
    Constructs a set A = {n <= N | n = 1 (mod 3)} and verifies
    that A+A contains no squares. It then prints the density of A.
    """
    print(f"Running the analysis for N = {N}")
    
    # Construct the set A
    A = [n for n in range(1, N + 1) if n % 3 == 1]

    # Pre-calculate squares up to the maximum possible sum 2*N
    max_sum = 2 * N
    squares = {k * k for k in range(1, int(math.sqrt(max_sum)) + 2)}

    # Verify that no sum of two elements in A is a square
    found_square_sum = False
    for i in range(len(A)):
        for j in range(i, len(A)):
            current_sum = A[i] + A[j]
            if current_sum in squares:
                print(f"Error: Found a square sum: {A[i]} + {A[j]} = {current_sum}")
                found_square_sum = True
                break
        if found_square_sum:
            break

    if not found_square_sum:
        print("Verification successful: A+A contains no squares for the constructed set A.")
        size_A = len(A)
        density = size_A / N
        print("The construction used is: A = {n in {1,...,N} | n = 1 (mod 3)}")
        print(f"The size of A is: |A| = {size_A}")
        print(f"The equation for the density c is: |A| / N = {size_A} / {N} = {density}")

# You can change the value of N to test for different sizes
N = 5000
solve_square_sum_problem(N)