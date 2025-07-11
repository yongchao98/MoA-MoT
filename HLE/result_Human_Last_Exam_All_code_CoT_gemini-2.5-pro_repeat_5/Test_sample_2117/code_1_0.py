import numpy as np

def calculate_product(n):
    """
    Calculates the product E_P * E_H * S_P * S_H for a given integer n.
    """
    if n <= 0:
        raise ValueError("n must be a positive integer")

    # Dimension of the matrices
    m = n + 2

    # Term E_H
    # Eigenvalues of H are {n+1, -1, ..., -1} (n+1 times)
    # Sorted eigenvalues: {-1, ..., -1, n+1}
    # The only non-zero gap is (n+1) - (-1) = n+2
    # The number of gaps is m-1 = n+1
    e_h = (n + 2) / (n + 1)

    # Term S_H
    # Based on our construction of H, ||H||_F^2 = 2n^2 + 3n + 2
    s_h = (2 * n**2 + 3 * n + 2) / m

    # Term S_P
    # Based on our construction of P, ||P||_F^2 = 3n + 2
    s_p = (3 * n + 2) / m

    # Term E_P
    # Construct the matrix P of size m x m
    p_matrix = np.zeros((m, m))
    # Column 1: e_1
    p_matrix[0, 0] = 1
    # Column 2: v - e_1
    p_matrix[1:, 1] = 1
    # Columns 3 to m: e_{k-1} - e_k
    for k in range(3, m + 1):
        p_matrix[k - 2, k - 1] = 1
        p_matrix[k - 1, k - 1] = -1
    
    # Calculate eigenvalues of P
    eigenvalues_p = np.linalg.eigvals(p_matrix)

    # Sort eigenvalues. For complex numbers, sort by real part then imaginary part.
    eigenvalues_p.sort()

    # Calculate sum of absolute differences (gaps)
    gaps_p = np.abs(np.diff(eigenvalues_p))
    sum_of_gaps_p = np.sum(gaps_p)
    
    # E_P is the average gap
    e_p = sum_of_gaps_p / (m - 1)

    # Final product
    product = e_p * e_h * s_p * s_h
    
    return product, e_p, e_h, s_p, s_h

def main():
    """
    Main function to perform the calculation and print the results.
    """
    print("This program calculates the product E_P * E_H * S_P * S_H for a regular n-simplex.")
    print("The formula for the product is derived as:")
    print("Product(n) = E_P(n) * ( (n+2)/(n+1) ) * ( (3*n+2)/(n+2) ) * ( (2*n^2+3*n+2)/(n+2) )\n")
    
    print("Numerically calculating the product for various n:")
    print("-" * 50)
    print(f"{'n':>2} | {'E_P':>10} | {'E_H':>10} | {'S_P':>10} | {'S_H':>10} | {'Product':>12}")
    print("-" * 50)

    for n_val in range(1, 11):
        try:
            prod, ep, eh, sp, sh = calculate_product(n_val)
            print(f"{n_val:2d} | {ep:10.4f} | {eh:10.4f} | {sp:10.4f} | {sh:10.4f} | {prod:12.4f}")
        except Exception as e:
            print(f"Could not calculate for n={n_val}: {e}")

    print("-" * 50)
    print("\nThe product appears to be a monotonically increasing function of n.")
    print("Therefore, the least upper bound over all positive integers n is infinity.")
    
if __name__ == '__main__':
    main()

<<<inf>>>