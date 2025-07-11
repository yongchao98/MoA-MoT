import numpy as np
from scipy.integrate import quad
import math

def calculate_primes_and_dimensions():
    """
    Calculates the dimensions d1 and d2 required for the integral.
    It does this by first generating a list of primes up to a sufficient limit.
    """
    # The largest prime index needed is 10231 for d2.
    # Estimate the size of the 10231st prime: n*ln(n) ~ 10231 * ln(10231) ~ 94440
    # Let's sieve up to 120,000 to be safe.
    limit = 120000
    primes = []
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for p in range(2, int(math.sqrt(limit)) + 1):
        if is_prime[p]:
            for i in range(p * p, limit + 1, p):
                is_prime[i] = False
    for p in range(2, limit + 1):
        if is_prime[p]:
            primes.append(p)

    # Get the required primes by their indices (1-based index)
    p_781 = primes[781 - 1]
    p_2321 = primes[2321 - 1]
    p_8231 = primes[8231 - 1]
    p_10231 = primes[10231 - 1]

    n1, p1 = p_8231, p_781
    n2, p2 = p_10231, p_2321

    # Calculate the dimensions
    # dim(M(n,p)) = n*p - p*(p+1)/2
    # Use integer arithmetic to maintain precision
    d1 = n1 * p1 - (p1 * (p1 + 1)) // 2
    d2 = n2 * p2 - (p2 * (p2 + 1)) // 2
    
    return d1, d2, n1, p1, n2, p2

def integrand(x, d1, d2):
    """
    The full integrand from the problem statement.
    """
    # Handle the removable singularity at x=0
    if x == 0:
        return 0.0

    # The second term of the integrand, which is x * exp(-x)
    integral_part_2 = x * math.exp(-x)

    # The first term of the integrand.
    # It involves very large powers, which must be handled carefully to avoid overflow.
    # For large d, 1/(1+x^(2d)) is effectively a step function (1 for x<1, 0 for x>1).
    # The difference of two such functions is 0 everywhere except for a tiny region
    # around x=1, which is numerically negligible.
    
    # We use logarithms to handle the large exponents without overflow.
    log_x = math.log(x)
    
    # Calculate h(x) = 1/(1+x^(2d2)) - 1/(1+x^(2d1))
    try:
        # Check for overflow before calling exp
        log_power_d1 = 2 * d1 * log_x
        if log_power_d1 > 700: # Threshold for float64 overflow in exp
            term_d1 = 0.0
        else:
            term_d1 = 1.0 / (1.0 + math.exp(log_power_d1))
    except OverflowError: # This case happens if 2*d1*log_x is huge
        term_d1 = 0.0

    try:
        log_power_d2 = 2 * d2 * log_x
        if log_power_d2 > 700:
            term_d2 = 0.0
        else:
            term_d2 = 1.0 / (1.0 + math.exp(log_power_d2))
    except OverflowError:
        term_d2 = 0.0
        
    h_x = term_d2 - term_d1
    
    # If h_x is effectively zero, we don't need to calculate the denominator
    if h_x == 0.0:
        integral_part_1 = 0.0
    else:
        # Use math.expm1 for better precision for x near 0
        sqrt_term = math.sqrt(math.expm1(2 * x))
        if sqrt_term == 0.0:
            integral_part_1 = 0.0
        else:
            integral_part_1 = h_x / (x * sqrt_term)

    return integral_part_1 + integral_part_2


def main():
    """
    Main function to solve the problem.
    """
    # Part 1: Calculate the sum S
    # l(n,p) is the injectivity radius of the Stiefel manifold M(n,p).
    # For n > p, this is known to be pi.
    # The prime indices ensure n > p for all terms in the sum.
    # The sum is over a 10x10 grid, so S = 10 * 10 * pi.
    S = 100 * np.pi

    # Part 2: Calculate the integral I
    # First, get the dimensions d1 and d2
    d1, d2, n1, p1, n2, p2 = calculate_primes_and_dimensions()

    # Numerically integrate the function.
    # The integral of the first part is expected to be 0, and the second part is 1.
    # So the total integral I should be 1.0.
    I, err = quad(integrand, 0, np.inf, args=(d1, d2))

    # Part 3: Calculate the final result
    result = S * I

    # Print the results as requested
    print("Step 1: Evaluating the sum")
    sum_str = "sum_{i=1}^{10} sum_{j=1}^{10} l(p_{(21367+i)}, p_{(14567+j)})"
    print(f"The sum is S = {sum_str} = 100 * pi")
    print(f"Value of S = {S:.10f}\n")
    
    print("Step 2: Evaluating the integral")
    print("The dimensions are:")
    print(f"d1 = dim(M(p_(8231), p_(781))) = dim(M({n1}, {p1})) = {d1}")
    print(f"d2 = dim(M(p_(10231), p_(2321))) = dim(M({n2}, {p2})) = {d2}")
    print(f"The integral evaluates to I = {I:.10f} (numerical error estimate: {err:.2e})\n")
    
    print("Step 3: Calculating the final result")
    print("The final result is S * I.")
    print("\nFinal equation with calculated values:")
    print(f"( {S} ) * ( {I} ) = {result}")
    
    print(f"\nFinal Answer: {result:.10f}")

if __name__ == "__main__":
    main()
<<<314.1592653589793>>>