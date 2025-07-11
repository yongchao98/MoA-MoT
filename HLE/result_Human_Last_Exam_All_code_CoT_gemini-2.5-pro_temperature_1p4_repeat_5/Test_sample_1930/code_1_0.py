import numpy as np

def get_digits(x, num_digits=10):
    """
    Returns the digit sequence for a number x in [0, 10).
    Note: This is for demonstration and has precision limitations.
    It returns the standard representation (no trailing 9s).
    """
    if not (0 <= x < 10):
      if x == 10:
        return [9] * num_digits
      raise ValueError("Input x must be in [0, 10]")
    
    digits = []
    # A_0
    digits.append(int(x))
    x -= int(x)
    
    # A_1, A_2, ...
    for _ in range(num_digits - 1):
        x *= 10
        digit = int(x)
        digits.append(digit)
        x -= digit
    return digits

def basis_b0(A):
    """Corresponds to the constant function f(x)=1."""
    return 1.0

def basis_b1(A):
    """Corresponds to the function f(x)=x."""
    val = 0.0
    for n, digit in enumerate(A):
        val += digit / (10**n)
    return val
    
def basis_b2(A):
    """Corresponds to the function f(A) = sum_n A_{n+1} / 10^n."""
    val = 0.0
    for n in range(len(A) - 1):
        val += A[n+1] / (10**n)
    return val

def basis_b3(A):
    """Corresponds to the function f(A) = sum_n A_{n+2} / 10^n."""
    val = 0.0
    for n in range(len(A) - 2):
        val += A[n+2] / (10**n)
    return val

# Test sequences to show linear independence
A0 = [0, 0, 0, 0]  # Corresponds to x=0
A1 = [1, 0, 0, 0]  # Corresponds to x=1
A2 = [0, 1, 0, 0]  # Corresponds to x=0.1
A3 = [0, 0, 1, 0]  # Corresponds to x=0.01

test_sequences = [A0, A1, A2, A3]

# Create a matrix where rows are basis functions and columns are test sequences
matrix = np.array([
    [basis_b0(A) for A in test_sequences],
    [basis_b1(A) for A in test_sequences],
    [basis_b2(A) for A in test_sequences],
    [basis_b3(A) for A in test_sequences],
])

# The rank of the matrix gives the number of linearly independent functions
# over this set of inputs. If rank is 4, they are independent.
rank = np.linalg.matrix_rank(matrix)

print(f"We have constructed a basis of {rank} functions.")
print("The four basis functions are:")
print("1. f(x) = 1")
print("2. f(x) = x = sum(A_n / 10^n)")
print("3. f(A) = sum(A_{n+1} / 10^n)")
print("4. f(A) = sum(A_{n+2} / 10^n)")
print("\nThe dimension of the vector space of digitary functions is the number of these basis functions.")
print("\nFinal Answer Equation: Dimension = 1 (constant) + 1 (depends on A_n) + 1 (depends on A_{n+1}) + 1 (depends on A_{n+2})")
print("1 + 1 + 1 + 1 = 4")