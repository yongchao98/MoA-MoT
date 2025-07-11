import math

def solve_circle_pdf():
    """
    Solves the problem by finding the circle equation, the PDF,
    the normalization constant alpha, and the probability P(X < 3).
    """
    # Step 1: Find the circle parameters
    # The perpendicular bisector of A(1,0) and B(10,0) is x = 5.5, so h=5.5.
    h = 5.5
    
    # Using point A(1,0): (1 - 5.5)^2 + (0 - k)^2 = r^2  => 20.25 + k^2 = r^2
    # Using point C(4,3): (4 - 5.5)^2 + (3 - k)^2 = r^2  => 2.25 + 9 - 6k + k^2 = r^2
    # Equating them: 20.25 + k^2 = 11.25 - 6k + k^2  => 9 = -6k => k = -1.5
    k = -1.5
    
    # r^2 = 20.25 + (-1.5)^2 = 20.25 + 2.25 = 22.5
    r2 = 22.5
    r = math.sqrt(r2)

    # Step 2: The function f(x) is the upper part of the circle y = k + sqrt(r^2 - (x-h)^2)
    # f(x) = -1.5 + sqrt(22.5 - (x-5.5)^2)

    # Helper function for the antiderivative of the circular part: integral of sqrt(r^2 - u^2)
    def G(x, h_val, r_val):
        u = x - h_val
        # The value inside arcsin must be in [-1, 1]
        arg_asin = u / r_val
        if arg_asin > 1.0: arg_asin = 1.0
        if arg_asin < -1.0: arg_asin = -1.0
        # The value inside sqrt must be non-negative
        arg_sqrt = r_val**2 - u**2
        if arg_sqrt < 0: arg_sqrt = 0
        
        term1 = (u / 2) * math.sqrt(arg_sqrt)
        term2 = (r_val**2 / 2) * math.asin(arg_asin)
        return term1 + term2

    # Step 3: Calculate the total integral to find alpha
    # Integral of the constant part: k * (10 - 1)
    integral_const_total = k * (10 - 1)
    # Integral of the circular part: G(10) - G(1)
    integral_circ_total = G(10, h, r) - G(1, h, r)
    # Total integral of f(x) from 1 to 10
    total_integral = integral_const_total + integral_circ_total
    
    # Calculate alpha
    alpha = 1 / total_integral

    # Step 4: Calculate the integral for P(X < 3)
    # The probability is defined from the start of the domain [1,10]
    # Integral of the constant part from 1 to 3: k * (3 - 1)
    integral_const_prob = k * (3 - 1)
    # Integral of the circular part from 1 to 3: G(3) - G(1)
    integral_circ_prob = G(3, h, r) - G(1, h, r)
    # Integral of f(x) from 1 to 3
    prob_integral = integral_const_prob + integral_circ_prob
    
    # Calculate the probability P(X < 3)
    prob_x_lt_3 = alpha * prob_integral
    
    # Outputting the numbers in the final equation as requested
    print("The final probability is calculated as the ratio of two integrals:")
    print(f"P(X < 3) = (Integral of f(x) from 1 to 3) / (Integral of f(x) from 1 to 10)")
    print(f"Value of Integral from 1 to 3: {prob_integral}")
    print(f"Value of Integral from 1 to 10: {total_integral}")
    print("-" * 30)
    print("Final Results:")
    print(f"The value of alpha is: {alpha}")
    print(f"The value of P(X < 3) is: {prob_x_lt_3}")


solve_circle_pdf()