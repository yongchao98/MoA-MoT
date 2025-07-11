import math

def solve_problem():
    """
    Solves the problem by finding alpha and P(X<3).
    """

    # Step 1: Define the circle's parameters.
    # From the points A=(1,0), B=(10,0), C=(4,3), we derive the circle's equation:
    # Center (h,k) and radius r.
    # (x-h)^2 + (y-k)^2 = r^2
    # The perpendicular bisector of AB is x=5.5, so h=5.5.
    # (1-5.5)^2 + (0-k)^2 = r^2 => 20.25 + k^2 = r^2
    # (4-5.5)^2 + (3-k)^2 = r^2 => 2.25 + 9 - 6k + k^2 = r^2
    # Solving 20.25 = 11.25 - 6k => 9 = -6k => k = -1.5
    # r^2 = 20.25 + (-1.5)^2 = 22.5
    # The function f(x) is the upper arc of the circle:
    # f(x) = k + sqrt(r^2 - (x-h)^2)
    h = 5.5
    k = -1.5
    r_squared = 22.5
    r = math.sqrt(r_squared)

    def anti_derivative_sqrt_part(u, r_val, r_sq_val):
        """
        Calculates the value of the anti-derivative of sqrt(r^2 - u^2).
        The anti-derivative is (u/2)*sqrt(r^2-u^2) + (r^2/2)*asin(u/r).
        """
        # Clamp the argument of asin to [-1, 1] to avoid domain errors from floating point inaccuracies
        asin_arg = max(-1.0, min(1.0, u / r_val))
        
        # Handle cases where u is very close to r or -r
        if abs(abs(u) - r_val) < 1e-9:
             sqrt_term = 0
        else:
             sqrt_term = (u / 2.0) * math.sqrt(r_sq_val - u**2)

        asin_term = (r_sq_val / 2.0) * math.asin(asin_arg)
        return sqrt_term + asin_term

    def calculate_integral_f(x_start, x_end):
        """
        Calculates the definite integral of f(x) from x_start to x_end.
        Integral of f(x) = Integral of (k + sqrt(r^2 - (x-h)^2)) dx
                         = Integral of k dx + Integral of sqrt(r^2 - (x-h)^2) dx
        """
        # First part: integral of constant k
        integral_k = k * (x_end - x_start)

        # Second part: integral of the sqrt term using substitution u = x-h
        u_start = x_start - h
        u_end = x_end - h
        
        integral_sqrt = anti_derivative_sqrt_part(u_end, r, r_squared) - \
                        anti_derivative_sqrt_part(u_start, r, r_squared)
                        
        return integral_k + integral_sqrt

    # Step 2: Calculate the total area under f(x) from 1 to 10 to find alpha.
    total_area = calculate_integral_f(1, 10)
    
    # alpha is the reciprocal of the total area.
    if total_area > 0:
        alpha = 1.0 / total_area
    else:
        alpha = float('inf') # Should not happen based on problem description

    # Step 3: Calculate the area under f(x) from 1 to 3 to find P(X<3).
    area_lt_3 = calculate_integral_f(1, 3)

    # P(X<3) is alpha times this area.
    prob_lt_3 = alpha * area_lt_3

    print(f"The value of alpha is: {alpha}")
    print(f"The value of P(X<3) is: {prob_lt_3}")

solve_problem()