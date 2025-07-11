import sympy

def solve_least_upper_bound():
    """
    This function calculates the least upper bound of the product E_P*E_H * S_P*S_H.

    The derivation is as follows:
    1.  Let n be a positive integer representing the dimension of the regular simplex.
    2.  The Cayley-Menger matrix C for a regular n-simplex with unit side length has size (n+1)x(n+1).
        Its eigenvalues are n (with multiplicity 1) and -1 (with multiplicity n).
    3.  We choose the eigendecomposition C = P*H*P^-1, where H is the diagonal matrix of eigenvalues and P is
        the orthogonal matrix of eigenvectors. This is a valid form of Hessenberg decomposition.
    4.  H is a diagonal matrix with entries {n, -1, ..., -1}.
    5.  The average eigenvalue gap of H, E_H = (max_eig - min_eig) / (size - 1)
       E_H = (n - (-1)) / ((n+1) - 1) = (n+1)/n.
    6.  The mean square of singular values of H, S_H = (1/size) * ||H||_F^2
       S_H = (1/(n+1)) * (n**2 + n*(-1)**2) = (n**2 + n)/(n+1) = n.
    7.  P is an orthogonal matrix. Its singular values are all 1.
    8.  The mean square of singular values of P, S_P = (1/size) * ||P||_F^2 = (1/(n+1)) * (n+1)*1**2 = 1.
    9.  The average eigenvalue gap of P, E_P = (max_eig_real_part - min_eig_real_part) / (size - 1).
        For any orthogonal matrix, this gap is at most 2. So, E_P <= 2/n.
    10. The total product is E_P * E_H * S_P * S_H = E_P * ((n+1)/n) * 1 * n = (n+1) * E_P.
    11. The product is thus <= (n+1) * (2/n) = 2*(n+1)/n.
    12. This function f(n) = 2*(n+1)/n is a decreasing function for n > 0.
    13. The least upper bound (supremum) is its value at n=1.
    14. At n=1, the value is 2*(1+1)/1 = 4. We showed in the reasoning that this value is achievable.
    """
    n = sympy.Symbol('n', positive=True, integer=True)

    print("Step 1: Define quantities as functions of n.")
    # E_H = (n+1)/n
    E_H = (n + 1) / n
    print(f"E_H = (n - (-1)) / ((n+1) - 1) = {E_H}")
    
    # S_H = n
    S_H = n
    print(f"S_H = (n^2 + n) / (n+1) = {S_H}")

    # S_P = 1
    S_P = 1
    print(f"S_P = 1 (since P is orthogonal)")
    
    # E_P <= 2/n
    # We will use the upper bound for E_P to find the upper bound for the product.
    E_P_bound = 2 / n
    print(f"E_P has an upper bound of 2/n")

    print("\nStep 2: Calculate the product's upper bound as a function of n.")
    product_bound = E_P_bound * E_H * S_P * S_H
    
    print(f"Product = E_P * E_H * S_P * S_H")
    print(f"Product <= (2/n) * ((n+1)/n) * 1 * n")
    
    # Simplify the expression for the product bound
    final_expression = sympy.simplify(product_bound)
    print(f"Simplified upper bound f(n) = {final_expression}")

    print("\nStep 3: Find the least upper bound of f(n) for n >= 1.")
    # The function f(n) = 2*(n+1)/n is decreasing for n > 0.
    # The least upper bound is the value at n=1.
    lub = sympy.limit(final_expression, n, 1)
    
    # We demonstrate that the bound is decreasing by checking values
    val_n1 = final_expression.subs(n, 1)
    val_n2 = final_expression.subs(n, 2)
    val_n3 = final_expression.subs(n, 3)
    
    print(f"f(1) = {val_n1}")
    print(f"f(2) = {val_n2}")
    print(f"f(3) = {val_n3}")
    
    print("\nAs f(n) is a decreasing function, its least upper bound is its value at n=1.")
    print(f"The least upper bound is {lub}.")
    
    # Final answer in requested format
    # The equation for n=1: E_P=2, E_H=2, S_P=1, S_H=1
    print("\nFor n=1, we showed that the bound is achieved:")
    print("E_P = (1 - (-1)) / 1 = 2")
    print("E_H = (1 - (-1)) / 1 = 2")
    print("S_P = 1")
    print("S_H = 1")
    print(f"Product at n=1: E_P*E_H*S_P*S_H = 2 * 2 * 1 * 1 = {val_n1}")

solve_least_upper_bound()