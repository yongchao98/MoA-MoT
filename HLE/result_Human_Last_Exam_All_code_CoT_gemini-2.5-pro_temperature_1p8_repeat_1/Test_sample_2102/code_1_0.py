import sympy as sp
import numpy as np

def get_taylor_coefficients(degree):
    """
    Computes the Taylor coefficients for the function (2/pi) * K(x) * exp(x)
    around x=0 up to a given degree.
    """
    x = sp.Symbol('x')
    # The complete elliptic integral of the first kind is K(m). The problem uses K(x),
    # so we use x as the parameter m.
    func = (2 / sp.pi) * sp.elliptic_k(x) * sp.exp(x)
    
    # sympy.series(f, x, x0, n) computes the series up to the (n-1)-th power.
    # To get coefficients up to degree `degree`, we need n = `degree` + 1.
    poly = func.series(x, 0, degree + 1)
    
    # Extract coefficients from the polynomial
    coeffs = [poly.coeff(x, i) for i in range(degree + 1)]
    
    # Convert from sympy's arbitrary-precision Float to standard Python float
    return [float(c) for c in coeffs]

# Main logic to find the smallest n
n = 1
while True:
    # For an n x n Hankel matrix S_n, where S_n[i,j] = c_{i+j}, we need
    # coefficients from c_0 up to c_{ (n-1) + (n-1) } = c_{2n-2}.
    # So the highest degree needed for the Taylor polynomial is 2n-2.
    required_degree = 2 * n - 2
    
    taylor_coeffs = get_taylor_coefficients(required_degree)
    
    # Construct the n x n Hankel matrix S_n.
    # The constructor np.linalg.hankel(c, r) takes the first column and the last row.
    # For our S_n:
    # First column is c = (c_0, c_1, ..., c_{n-1})
    # Last row is r = (c_{n-1}, c_n, ..., c_{2n-2})
    if n == 1:
        # For n=1, S_1 is simply the 1x1 matrix [[c_0]]
        S_n = np.array([[taylor_coeffs[0]]])
    else:
        first_col = taylor_coeffs[0:n]
        last_row = taylor_coeffs[n-1 : 2*n-1]
        S_n = np.linalg.hankel(first_col, r=last_row)
        
    # The Hankel matrix from real coefficients is real and symmetric.
    # The eigenvalues of W_n are the same as the eigenvalues of S_n.
    # We use np.linalg.eigvalsh for efficient eigenvalue computation for Hermitian matrices.
    eigenvalues = np.linalg.eigvalsh(S_n)
    
    # f(n) is the sum of the absolute cubes of the eigenvalues.
    f_n = np.sum(np.abs(eigenvalues)**3)
    
    # Check if the condition f(n) > 10 is met.
    if f_n > 10:
        # We found the smallest n. Now calculate the final result.
        
        # For a symmetric matrix S_n, its Weyr form W_n is a diagonal matrix of its eigenvalues.
        # The infinity norm ||W_n||_inf is the maximum of the absolute values of the
        # diagonal elements (the eigenvalues).
        infinity_norm_W_n = np.max(np.abs(eigenvalues))
        
        # Calculate the final product.
        result = n * infinity_norm_W_n
        
        # Print the required output equation and exit the loop.
        print(f"{n} * {infinity_norm_W_n} = {result}")
        
        # Return the final answer in the specified format
        print(f'<<<{result}>>>')

        break
        
    # If the condition is not met, increment n and continue the search.
    n += 1