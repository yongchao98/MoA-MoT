import numpy as np
from scipy.integrate import quad

def sinc(x):
    """
    Defines the unnormalized sinc function, sin(x)/x.
    """
    if x == 0:
        return 1.0
    return np.sin(x) / x

def integrand_factory(n):
    """
    Returns the integrand function for a given n.
    The integrand is the product of sinc functions.
    """
    def integrand(x):
        prod = 1.0
        for k in range(1, n + 1):
            # Using x/k can cause division by zero if x is 0 at k=0, but our k starts at 1.
            prod *= sinc(x / k)
        return prod
    return integrand

# The value that I_n is compared to.
pi_2 = np.pi / 2

print(f"Analyzing the Borwein Integral I_n = ∫[0,∞] Π[k=1 to n] sinc(x/k) dx")
print(f"The proposition P(n) is I_n = π/2, where π/2 ≈ {pi_2:.12f}")
print("-" * 65)
print("n |      I_n      |    I_n - π/2    |  Condition: Σ_{k=2 to n} 1/k")
print("-" * 65)

# This sum determines when the integral stops being π/2. The condition is Σ <= 1.
sum_of_reciprocals = 0.0

# Calculate and print values for n=1 to 8.
for n in range(1, 9):
    # Get the integrand for the current n
    integrand = integrand_factory(n)
    
    # Perform numerical integration from 0 to infinity
    # The 'limit' argument increases subdivisions for better accuracy if needed
    val, err = quad(integrand, 0, np.inf, limit=200)
    
    # Update the sum for the condition check
    if n > 1:
        sum_of_reciprocals += 1.0/n
        
    # Print the results in a formatted table
    # The final column shows the sum. The condition for I_n = π/2 is that this sum is <= 1.
    print(f"{n:1d} | {val:15.12f} | {val - pi_2:15.8e} | {sum_of_reciprocals:15.4f}")
