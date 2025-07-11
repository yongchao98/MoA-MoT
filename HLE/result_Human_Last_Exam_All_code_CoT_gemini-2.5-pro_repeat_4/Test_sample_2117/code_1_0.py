import numpy as np

def solve():
    """
    Calculates the least upper bound of the product E_P*E_H * S_P*S_H.

    The overall product simplifies to the function f(n) = 2*(n+2)/(n+1).
    This function is maximized for the smallest positive integer n, which is n=1.
    This script calculates the components of the product for n=1 and prints the final result.
    """
    n = 1
    m = n + 2

    # E_H: Average eigenvalue gap of H
    # Eigenvalues of C_n (and H) are (n+1) and -1 (with multiplicity n+1)
    # E_H = (lambda_max - lambda_min) / (m-1)
    lambda_max_H = n + 1
    lambda_min_H = -1
    E_H = (lambda_max_H - lambda_min_H) / (m - 1)

    # S_H: Mean square of singular values of H
    # For symmetric H, sigma_i = |lambda_i|. S_H = (1/m) * sum(|lambda_i|^2)
    # Sum of squares = (n+1)*|-1|^2 + |n+1|^2 = (n+1) + (n+1)^2 = (n+1)*(n+2)
    S_H = (n + 1) * (n + 2) / m

    # P is an orthogonal matrix from the standard decomposition of a symmetric matrix.
    # S_P: Mean square of singular values of P
    # For orthogonal P, all singular values are 1. So S_P = 1.
    S_P = 1.0

    # E_P: Average eigenvalue gap of P
    # For the given C_n, P has eigenvalues 1 (n+1 times) and -1 (once).
    # E_P = (lambda_max - lambda_min) / (m-1)
    lambda_max_P = 1
    lambda_min_P = -1
    E_P = (lambda_max_P - lambda_min_P) / (m - 1)

    # The final product
    product = E_P * E_H * S_P * S_H

    print(f"For n = {n}:")
    print(f"E_P = ({lambda_max_P} - ({lambda_min_P})) / ({m} - 1) = {E_P}")
    print(f"E_H = ({lambda_max_H} - ({lambda_min_H})) / ({m} - 1) = {E_H}")
    print(f"S_P = {S_P}")
    print(f"S_H = ({n+1}*1^2 + {n+1}^2) / {m} = {S_H}")
    print(f"Product = {E_P} * {E_H} * {S_P} * {S_H} = {product}")
    
solve()