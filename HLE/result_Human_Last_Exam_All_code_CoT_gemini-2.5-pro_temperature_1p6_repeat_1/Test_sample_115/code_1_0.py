import numpy as np
from scipy.linalg import hadamard

def get_A(n):
    """
    Generates the matrix A_n.
    The rows and columns are indexed by integers from 0 to 2^n - 1,
    which correspond to subsets of an n-element set.
    """
    size = 2**n
    A = np.zeros((size, size))
    subsets = []
    for i in range(size):
        s = set()
        for j in range(n):
            if (i >> j) & 1:
                s.add(j)
        subsets.append(s)
    
    for i in range(size):
        for j in range(size):
            if subsets[i].intersection(subsets[j]):
                A[i, j] = 1
    return A

def get_W(n):
    """
    Generates the Walsh-Hadamard matrix W_n.
    This is a normalized real unitary matrix.
    """
    size = 2**n
    W = hadamard(size) / np.sqrt(size)
    return W

def calculate_alpha():
    """
    Calculates the norms for n=1 and n=2 and their ratio.
    """
    print("Step 1: Calculate the spectral norm for n=1.")
    
    # n=1
    n1 = 1
    A1 = get_A(n1)
    W1 = get_W(n1)
    M1 = A1 * W1 # Hadamard product
    norm1 = np.linalg.norm(M1, ord=2)
    print("A_1 = \n", A1)
    print("W_1 = \n", W1)
    print("M_1 = A_1 * W_1 = \n", M1)
    print(f"The spectral norm for n=1 is ||A_1 o W_1|| = {norm1:.4f}")
    # Exact value is 1/sqrt(2)
    
    print("\nStep 2: Calculate the spectral norm for n=2.")
    # n=2
    n2 = 2
    A2 = get_A(n2)
    W2 = get_W(n2)
    M2 = A2 * W2 # Hadamard product
    norm2 = np.linalg.norm(M2, ord=2)
    print("A_2 = \n", A2)
    print("W_2 = \n", W2)
    print("M_2 = A_2 * W_2 = \n", M2)
    print(f"The spectral norm for n=2 is ||A_2 o W_2|| = {norm2:.4f}")
    # Exact value is sqrt(3)/2
    
    print("\nStep 3: Estimate alpha as the ratio of the norms.")
    alpha = norm2 / norm1
    print(f"The ratio of the norms is {norm2:.4f} / {norm1:.4f} = {alpha:.4f}")
    print("\nBased on these calculations, the value of alpha is conjectured to be sqrt(3/2).")
    print(f"sqrt(3/2) is approximately {np.sqrt(3/2):.4f}")
    print("\nThe equation for the growth rate of c_n is c_n = Theta((sqrt(3/2))^n).")
    print("\nFinal Answer:")
    # Using 'print' to show the numbers in the final expression
    numerator = 3
    denominator = 2
    exponent = "n"
    print(f"(sqrt({numerator}/{denominator}))^{exponent}")

calculate_alpha()
