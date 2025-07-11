import numpy as np

def solve():
    """
    This function calculates and verifies the deadbeat observer gain matrix L.
    """
    # Define the system matrices
    A = np.array([
        [-1, 0, 0,  1],
        [ 1, 0, 0,  2],
        [ 0, 1, 0, -1],
        [-1, 0, 1, -1]
    ])

    C = np.array([
        [1, 0, 0, 0],
        [1, 0, 0, 1]
    ])

    # The observer gain L is designed to make (A - LC)^3 = 0.
    # Through algebraic construction based on the structure of A and C,
    # we can find a suitable L.
    # We choose l31=0 and l41=0 for simplicity, as they are free parameters in this construction.
    L = np.array([
        [-2.0,  1.0],
        [-1.0,  2.0],
        [ 0.0, -1.0],
        [ 0.0, -1.0]
    ])

    # Verify the design
    # 1. Compute the error dynamics matrix
    A_o = A - L @ C

    # 2. The observability index is 3, so we check if (A - LC)^3 is the zero matrix.
    A_o_3 = np.linalg.matrix_power(A_o, 3)

    print("The system observability index is 3.")
    print("We need to design L such that (A - LC)^3 = 0.\n")
    print("The constructed observer gain matrix L is:")
    print(L)
    print("\nTo verify, we compute (A - LC)^3:")
    print(A_o_3)
    
    # The small non-zero values are due to floating-point arithmetic.
    is_zero = np.allclose(A_o_3, np.zeros_like(A_o_3))
    print(f"\nIs (A - LC)^3 a zero matrix? {is_zero}")


solve()

<<<[[ -2.   1.]
 [ -1.   2.]
 [  0.  -1.]
 [  0.  -1.]]>>>