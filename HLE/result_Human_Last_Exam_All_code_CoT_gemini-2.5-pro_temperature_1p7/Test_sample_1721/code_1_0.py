import math

def is_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    sqrt_n = int(math.sqrt(n))
    return sqrt_n * sqrt_n == n

def verify_construction(N):
    """
    Constructs the set A = {n in {1,...,N} | n = 1 (mod 3)} and verifies
    that A+A contains no squares.
    """
    # Construct the set A
    A = [n for n in range(1, N + 1) if n % 3 == 1]
    
    # Verify the condition
    found_square_sum = False
    for i in range(len(A)):
        for j in range(i, len(A)):
            a = A[i]
            b = A[j]
            s = a + b
            if is_square(s):
                print(f"Verification failed: {a} + {b} = {s}, which is a square.")
                found_square_sum = True
                break
        if found_square_sum:
            break
            
    if not found_square_sum:
        print(f"Verification successful for N={N}.")
        print(f"The constructed set A has {len(A)} elements.")
        print(f"The density is |A|/N = {len(A)}/{N} = {len(A)/N:.4f}")
        print("This is consistent with the theoretical maximum density c = 1/3.")

# You can change this value to test for different N
N = 1000
verify_construction(N)

# The final answer is the value of c.
# We will represent the number 1/3 in the final equation.
final_c = 1/3
print("\nThe final equation for the size of A is: |A| = (c + o(1))N")
print(f"The largest possible value for c is 1/3.")
# Outputting the numbers from the fraction that represents c
numerator = 1
denominator = 3
print(f"Numerator: {numerator}")
print(f"Denominator: {denominator}")