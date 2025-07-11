import numpy as np
from scipy.integrate import quad
from scipy.special import sinc as scipy_sinc # sinc(x) = sin(pi*x)/(pi*x)

# Define the integrand for I_n
def integrand(x, n):
    """
    Computes the product Π_{k=1 to n} sinc(x/k).
    The standard sinc(y) is sin(y)/y.
    The scipy.special.sinc(y) is defined as sin(πy)/(πy).
    To get sin(z)/z using scipy's sinc, we can use the relation:
    sin(z)/z = sin(π * (z/π)) / (π * (z/π)) = scipy_sinc(z/π).
    The argument to our sinc function in the integral is x/k.
    So we compute scipy_sinc( (x/k) / π ).
    """
    prod = 1.0
    for k in range(1, n + 1):
        prod *= scipy_sinc(x / (k * np.pi))
    return prod

# --- Main execution ---
# Header for the output
print("Numerical Evaluation of Borwein Integrals I_n")
print("-" * 60)
print(f"The reference value for the proposition P(n): I_n = π/2")
print(f"π/2 ≈ {np.pi/2:.10f}")
print("-" * 60)

# Iterate from n=1 to n=9 to observe the behavior of the integral
for n in range(1, 10):
    # Perform the numerical integration from 0 to infinity.
    # We increase the number of subdivisions for better accuracy on oscillatory integrals.
    integral_value, error = quad(integrand, 0, np.inf, args=(n,), limit=250)

    # The proposition P(n) states that I_n = π/2. Let's check it.
    # The 'equation' is I_n = π/2. The numbers are I_n's value and π/2's value.
    diff = integral_value - np.pi/2
    
    # Determine the relationship using a small tolerance for floating point comparisons
    if abs(diff) < 1e-9:
        relation_str = "is equal to"
        proposition_status = "True"
    elif diff < 0:
        relation_str = "is less than"
        proposition_status = "False"
    else:
        relation_str = "is greater than"
        proposition_status = "False"

    # Print the result for each n
    print(f"For n = {n}:")
    print(f"  The proposition is I_{n} = π/2.")
    print(f"  Numerically, I_{n} ≈ {integral_value:.10f}")
    print(f"  This value {relation_str} π/2.")
    print(f"  Therefore, P({n}) is numerically {proposition_status}.")
    print(f"  Difference (I_n - π/2) is {diff:+.3e}\n")
