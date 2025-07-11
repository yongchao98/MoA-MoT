import math

def solve_square_sum_problem(N):
    """
    This function demonstrates that the set A = {n | n=1 (mod 3)}
    has sums a+b that are never square numbers.
    It calculates the density of this set for a given N.
    """
    # Step 1: Construct the set A = {n in {1,...,N} | n = 1 (mod 3)}
    A = [n for n in range(1, N + 1) if n % 3 == 1]

    # Step 2: Verify that A+A contains no squares
    # For efficiency, pre-compute all squares up to the maximum possible sum 2*N
    max_sum = 2 * N
    squares = {k * k for k in range(1, int(math.sqrt(max_sum)) + 2)}

    is_violated = False
    for i in range(len(A)):
        for j in range(i, len(A)):
            a1 = A[i]
            a2 = A[j]
            s = a1 + a2
            if s in squares:
                print(f"Violation found: {a1} + {a2} = {s}, which is a square.")
                is_violated = True
                break
        if is_violated:
            break

    # Step 3: Print the results
    if not is_violated:
        print(f"Verification successful for N = {N}.")
        print("The set A = {n | 1 <= n <= N, n = 1 (mod 3)} was constructed.")
        print("The sumset A+A contains no square numbers.")
    else:
        print("Verification failed.")

    print(f"\nSize of set A: |A| = {len(A)}")
    print(f"Ratio |A|/N = {len(A) / N:.5f}")

    # Step 4: Output the final answer for c
    numerator = 1
    denominator = 3
    c = numerator / denominator

    print(f"\nThe largest number c is the fraction {numerator}/{denominator}.")
    print(f"The numbers in this fraction are {numerator} and {denominator}.")
    print(f"As a decimal, c = {c:.5f}...")


# Set a value for N to run the demonstration
N = 2000
solve_square_sum_problem(N)
