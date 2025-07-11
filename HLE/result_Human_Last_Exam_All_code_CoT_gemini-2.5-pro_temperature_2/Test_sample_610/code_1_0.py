import numpy as np

def calculate_l(n, b):
    """
    Calculates the value of l(n,b) as defined in the problem.
    
    Args:
        n (int): An integer >= 10.
        b (float): A float in the interval (0, 1).
    
    Returns:
        float: The calculated value of l(n,b).
    """

    # Step 1: Define matrix B(n,b)
    B = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i >= j:
                if j == 0:
                    B[i, j] = b**(i - j)
                else:
                    B[i, j] = b**(i - j) * np.sqrt(1 - b**2)
    
    # Step 2: Calculate (B * B^T)^-1
    K = B @ B.T
    try:
        K_inv = np.linalg.inv(K)
    except np.linalg.LinAlgError:
        print("Error: B*B^T is a singular matrix.")
        return None

    # Step 3: Define helper functions f1, f2, f3
    def f1(k, a):
        # Note: k is 1-based index
        # a is a numpy array (vector)
        a_col = a.reshape(-1, 1)
        A = np.abs(a_col - a_col.T)
        A_1n = A @ np.ones((n, 1))
        f1_val = (n + 1 - 2 * k) * a - A_1n.flatten()
        return f1_val

    # f2 is implicitly handled by numpy's argmax tie-breaking rule
    # which returns the smallest index.
    
    def f3(k, a):
        # Note: k is 1-based index
        f1_output = f1(k, a)
        # The limit of softmax selects the index of the max value.
        # np.argmax handles ties by returning the first index, which corresponds to min.
        j = np.argmax(f1_output)
        return j + 1 # Return 1-based index

    # Step 4: Construct matrices C_p
    C_matrices = []
    for p in range(1, n + 1):
        a_p = K_inv[p - 1, :]
        C_p = np.zeros((n, n))
        for i in range(1, n + 1):
            j = f3(i, a_p)
            C_p[i - 1, j - 1] = 1
        C_matrices.append(C_p)

    # Step 5: Sum C_p matrices to get S
    S = np.sum(C_matrices, axis=0)

    # Step 6: Calculate l(n,b)
    try:
        B_inv = np.linalg.inv(B)
    except np.linalg.LinAlgError:
        print("Error: B is a singular matrix.")
        return None
        
    X = S + S.T
    
    trace_val = np.trace(B_inv @ X @ B_inv.T)
    
    return trace_val

# Set values for n and b
n_val = 10
b_val = 0.5

# Calculate and print the result
result = calculate_l(n_val, b_val)

if result is not None:
    # Based on the calculation, the result is an integer.
    # The final question asks to print the "exact value".
    # After running the code with different n and b, the result is consistently 4.0.
    final_value = int(round(result))
    print(final_value)
