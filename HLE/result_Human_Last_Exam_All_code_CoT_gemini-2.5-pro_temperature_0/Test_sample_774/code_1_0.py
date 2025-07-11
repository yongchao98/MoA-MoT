import numpy as np
from scipy.integrate import quad
from sympy import prime

def get_dim(n_idx, p_idx):
    """Calculates the dimension of the manifold M(n,p)."""
    n = prime(n_idx)
    p = prime(p_idx)
    # The dimension formula is np - p(p+1)/2
    dim = n * p - (p * (p + 1)) // 2
    return n, p, dim

def calculate_sum_part():
    """Calculates the sum part of the expression."""
    # The injectivity radius l(n,p) of the Stiefel manifold is pi.
    # The sum is over a 10x10 grid.
    return 10 * 10 * np.pi

def integrand_I1(x, d1, d2):
    """
    Integrand for the first part of the integral (I1), handled carefully for numerical stability.
    I1 = integral of K(x) * (phi(x,d2) - phi(x,d1))
    where K(x) = 1/(x*sqrt(exp(2x)-1)) and phi(x,d) = 1/(1+x^2d)
    """
    # The kernel of the integrand
    # Handle singularity at x=0, where sqrt(e^2x - 1) ~ sqrt(2x)
    # The whole term is ~ x^(2d-3/2) which is integrable.
    if x == 0:
        return 0.0
    kernel = 1.0 / (x * np.sqrt(np.expm1(2 * x)))

    # The difference term, phi(x,d2) - phi(x,d1)
    # This term is numerically challenging due to the large exponents d1 and d2.
    # We split the computation based on x.
    if x == 1.0:
        diff = 0.0
    elif x < 1.0:
        # For x < 1, x^(2d) is safe to compute as it tends to 0.
        phi1 = 1.0 / (1.0 + x**(2 * d1))
        phi2 = 1.0 / (1.0 + x**(2 * d2))
        diff = phi2 - phi1
    else: # x > 1.0
        # For x > 1, x^(2d) would overflow. We use logarithms.
        # 1/(1+x^A) = 1/(1+exp(A*log(x))) = exp(-A*log(x))/(exp(-A*log(x))+1)
        # For large A, this is approximately exp(-A*log(x)).
        log_x = np.log(x)
        # The approximation is accurate enough as exp(-2*d*log_x) is extremely small.
        phi1_approx = np.exp(-2 * d1 * log_x)
        phi2_approx = np.exp(-2 * d2 * log_x)
        diff = phi2_approx - phi1_approx
        
    return kernel * diff

def calculate_integral_part(d1, d2):
    """Calculates the integral part of the expression."""
    # The integral I can be split into I = I1 + I2
    # I2 is integral of x*exp(-x) from 0 to inf, which is 1.
    I2 = 1.0
    
    # I1 is the complex integral which we expect to be 0.
    # We compute it numerically to verify.
    # We split the integral at x=1 for numerical stability.
    I1_part1, err1 = quad(integrand_I1, 0, 1, args=(d1, d2))
    I1_part2, err2 = quad(integrand_I1, 1, np.inf, args=(d1, d2))
    I1 = I1_part1 + I1_part2
    
    # The total integral value
    I = I1 + I2
    return I, I1, I2

# --- Main Calculation ---

# 1. Calculate dimensions d1 and d2
n1_idx, p1_idx = 8231, 781
n2_idx, p2_idx = 10231, 2321

n1, p1, d1 = get_dim(n1_idx, p1_idx)
n2, p2, d2 = get_dim(n2_idx, p2_idx)

# 2. Calculate the sum part
sum_val = calculate_sum_part()

# 3. Calculate the integral part
integral_val, I1_val, I2_val = calculate_integral_part(d1, d2)

# 4. Calculate the final result
final_result = sum_val * integral_val

# --- Output Results ---
print("--- Step-by-step Calculation ---")
print("\n1. Analysis of the Sum Part:")
print(f"The injectivity radius l(n,p) is pi. The sum is 10x10 terms.")
print(f"Value of the sum = 100 * pi = {sum_val}")

print("\n2. Analysis of the Integral Part:")
print("The integral is split into two parts: I = I1 + I2.")
print(f"The first dimension is d1 = dim(M(p_{n1_idx}, p_{p1_idx})) = dim(M({n1}, {p1})) = {d1}")
print(f"The second dimension is d2 = dim(M(p_{n2_idx}, p_{p2_idx})) = dim(M({n2}, {p2})) = {d2}")
print(f"The integral I2 = integral(x*exp(-x)) = {I2_val}")
print(f"The integral I1 is computed numerically and is expected to be 0.")
print(f"Numerical value of I1 = {I1_val:.6g}")
print(f"Total value of the integral I = I1 + I2 = {integral_val:.6g}")

print("\n3. Final Result:")
print("The final result is the product of the sum and the integral.")
print(f"Final Result = ({sum_val}) * ({integral_val})")
print(f"             = {final_result}")

print("\n--- Final Equation with Calculated Values ---")
print(f"( {sum_val} ) * ( {integral_val} ) = {final_result}")

<<<314.1592653589793>>>