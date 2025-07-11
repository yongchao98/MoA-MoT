import math
import sympy

def solve_problem():
    """
    Solves the given mathematical problem step-by-step.
    """
    # Part 1: The Summation
    # The term l(n,p) is the global injectivity radius of the Stiefel manifold M(n,p).
    # For the Stiefel manifold with the standard Euclidean metric, the injectivity radius is pi.
    # This is a standard result from differential geometry.
    # Therefore, l(n,p) = pi for all n, p >= 5.
    injectivity_radius = math.pi
    
    # The sum is a double summation over i from 1 to 10 and j from 1 to 10.
    # sum = sum_{i=1 to 10} sum_{j=1 to 10} l(p_(21367+i), p_(14567+i))
    # Since l(n,p) is a constant pi, the sum is 10 * 10 * pi.
    sum_val = 10 * 10 * injectivity_radius
    
    print("Step 1: Calculating the sum")
    print(f"The injectivity radius l(n,p) is a constant: pi")
    print(f"The sum is a 10x10 sum of pi, which is 100 * pi.")
    print(f"Value of the sum = {sum_val}\n")

    # Part 2: The Integral
    # The integral can be split into two parts:
    # Integral = Integral_A + Integral_B
    # Integral_A = integral from 0 to inf of (x^(2*d1) - x^(2*d2)) / (x * (1+x^(2*d2)) * (1+x^(2*d1)) * sqrt(e^(2x)-1)) dx
    # Integral_B = integral from 0 to inf of (x^2 * e^(-x) * (1+x^(2*d2))*(1+x^(2*d1))*sqrt(e^(2x)-1)) / (x * (1+x^(2*d2))*(1+x^(2*d1))*sqrt(e^(2x)-1)) dx

    # Let's analyze Integral_B first.
    # The integrand simplifies to x * e^(-x).
    # The integral of x*e^(-x) from 0 to inf is Gamma(2) = 1! = 1.
    integral_B_val = 1.0

    # Now, let's analyze Integral_A.
    # We need the dimensions d1 and d2.
    # dim(M(n,p)) = n*p - p*(p+1)/2
    
    # Find the required prime numbers
    n1_idx, p1_idx = 8231, 781
    n2_idx, p2_idx = 10231, 2321
    
    n1 = sympy.prime(n1_idx)
    p1 = sympy.prime(p1_idx)
    n2 = sympy.prime(n2_idx)
    p2 = sympy.prime(p2_idx)

    # Define the dimension function
    def dim_manifold(n, p):
        return n * p - p * (p + 1) // 2

    d1 = dim_manifold(n1, p1)
    d2 = dim_manifold(n2, p2)

    print("Step 2: Calculating the integral")
    print("The integral is composed of two parts, Integral_A and Integral_B.")
    print(f"The second part, Integral_B, simplifies to the integral of x*exp(-x) from 0 to infinity, which is {integral_B_val}.")
    
    print("\nThe first part, Integral_A, depends on manifold dimensions d1 and d2.")
    print(f"n1 = p_({n1_idx}) = {n1}, p1 = p_({p1_idx}) = {p1}")
    print(f"d1 = dim(M(n1, p1)) = {d1}")
    print(f"n2 = p_({n2_idx}) = {n2}, p2 = p_({p2_idx}) = {p2}")
    print(f"d2 = dim(M(n2, p2)) = {d2}")

    # The exponents d1 and d2 are enormous, making numerical integration impossible.
    # The term (x^(2*d1) - x^(2*d2)) / ((1+x^(2*d1))*(1+x^(2*d2))) can be rewritten as
    # 1/(1+x^(2*d2)) - 1/(1+x^(2*d1)).
    # For such large d, the function 1/(1+x^(2d)) is extremely close to a step function that drops from 1 to 0 at x=1.
    # The difference between two such step-like functions (for d1 and d2) is a function that is non-zero only in an infinitesimally small region around x=1.
    # Therefore, its integral against a smooth function is effectively zero.
    integral_A_val = 0.0
    
    print(f"\nAs d1 and d2 are extremely large, the integrand of Integral_A is effectively zero everywhere.")
    print(f"Thus, the value of Integral_A is taken to be {integral_A_val}.")

    integral_val = integral_A_val + integral_B_val
    print(f"Total value of the integral = {integral_A_val} + {integral_B_val} = {integral_val}\n")

    # Part 3: Final Calculation
    final_result = sum_val * integral_val
    
    print("Step 3: Final Calculation")
    print(f"Final Result = (Value of Sum) * (Value of Integral)")
    print(f"Final Result = {sum_val} * {integral_val} = {final_result}")

solve_problem()
<<<314.1592653589793>>>