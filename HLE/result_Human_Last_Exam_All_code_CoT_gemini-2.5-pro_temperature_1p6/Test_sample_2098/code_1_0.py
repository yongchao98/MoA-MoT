import numpy as np
from scipy.integrate import quad

def solve_problem():
    """
    Solves the problem by finding the coefficients for y1(x), determining n,
    and calculating the final integral.
    """
    # Step 1: Find the coefficients for y1(x) = C1/x^3 + C2/x^4 + C3/x^5
    A = np.array([
        [1/2**3, 1/2**4, 1/2**5],
        [1/6**3, 1/6**4, 1/6**5],
        [1/10**3, 1/10**4, 1/10**5]
    ])
    b = np.array([667, 2/9, 1/625])
    
    coeffs = np.linalg.solve(A, b)
    C1, C2, C3 = coeffs
    
    print(f"The coefficients for y1(x) are:")
    print(f"C1 = {C1:.0f}")
    print(f"C2 = {C2:.0f}")
    print(f"C3 = {C3:.0f}")
    print("-" * 20)

    # The function y1(x)
    def y1(x):
        return C1/x**3 + C2/x**4 + C3/x**5

    # Step 2: Determine the minimal integer 'n' for non-intersection
    # The intersection condition is y1(x)/x = 1/n
    # Let g(x) = y1(x)/x
    def g(x):
        return C1/x**4 + C2/x**5 + C3/x**6

    # We need to find the range of g(x) for x > 0.
    # Find critical points by finding roots of g'(x)=0.
    # g'(x) = -4*C1/x^5 - 5*C2/x^6 - 6*C3/x^7
    # g'(x) = 0 => -4*C1*x^2 - 5*C2*x - 6*C3 = 0
    # A quadratic equation for x.
    p = [-4 * C1, -5 * C2, -6 * C3]
    roots = np.roots(p)
    
    # We are interested in positive real roots
    positive_real_roots = [r.real for r in roots if r.imag == 0 and r.real > 0]
    
    # g(x) -> +inf as x->0+, g(x) -> 0 as x->inf.
    # Let's find the minimum value of g(x).
    min_g = float('inf')
    critical_x = -1
    for r in positive_real_roots:
        val = g(r)
        if val < min_g:
            min_g = val
            critical_x = r

    # For non-intersection, 1/n must be outside the range [min_g, inf)
    # So we need 1/n < min_g
    # This implies n > 1/min_g
    
    n_float = 1 / min_g
    minimal_n = int(np.floor(n_float)) + 1
    
    print(f"The critical point of g(x) = y1(x)/x is at x = {critical_x:.4f}")
    print(f"The minimum value of g(x) is {min_g:.4f}")
    print(f"Condition for non-intersection: n > {n_float:.4f}")
    print(f"The minimal integer 'n' is {minimal_n}")
    print("-" * 20)

    # Step 3 & 4: Calculate the integral
    # As argued in the plan, the integration interval is [2, 10].
    
    # Define the indefinite integral of y1(x)
    # integral(C1*x^-3 + C2*x^-4 + C3*x^-5) dx = C1*x^-2/-2 + C2*x^-3/-3 + C3*x^-4/-4
    def integral_y1(x):
        return C1*x**-2/-2 + C2*x**-3/-3 + C3*x**-4/-4

    integral_val_at_10 = integral_y1(10)
    integral_val_at_2 = integral_y1(2)
    
    result = integral_val_at_10 - integral_val_at_2

    # Verification using quad
    # result_quad, _ = quad(y1, 2, 10)
    
    print("The final calculation for the integral:")
    # Using fractions for precision in display
    # C1=564, C2=-9416, C3=37920
    # Antiderivative: -282/x^2 + 9416/(3*x^3) - 9480/x^4
    # At x=10: -282/100 + 9416/3000 - 9480/10000 = -2.82 + 3.13866... - 0.948 = -0.629333
    # At x=2: -282/4 + 9416/24 - 9480/16 = -70.5 + 392.333... - 592.5 = -270.666...
    # Result = (-0.629333...) - (-270.666...) = 270.037333...
    
    print(f"Integral value at x=10: {-282}/{10**2} + {9416}/(3*{10**3}) - {9480}/{10**4} = {integral_val_at_10:.4f}")
    print(f"Integral value at x=2: {-282}/{2**2} + {9416}/(3*{2**3}) - {9480}/{2**4} = {integral_val_at_2:.4f}")
    print(f"Result = ({integral_val_at_10:.4f}) - ({integral_val_at_2:.4f}) = {result:.4f}")
    
    final_answer = 101264/375
    print(f"The exact result is 101264/375, which is approximately {final_answer:.4f}")
    print(f"<<<{final_answer}>>>")


solve_problem()