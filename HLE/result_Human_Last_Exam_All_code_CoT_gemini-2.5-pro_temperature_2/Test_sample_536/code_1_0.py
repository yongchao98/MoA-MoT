import math

def solve_for_inner_product():
    """
    Calculates the limit of the inner product <b_p, z_p> based on the given geometric relations.
    """
    
    # Given values
    # <h_p, b_p>
    c_hb = 0.9375
    # <h_p, z_p>
    c_hz = 0.9
    
    # We want to find x = <b_p, z_p>.
    # The derived quadratic equation is: x^2 - (2*c_hb*c_hz)*x + (c_hb^2 + c_hz^2 - 1) = 0
    # Coefficients for ax^2 + bx + c = 0
    a = 1
    b = -2 * c_hb * c_hz
    c = c_hb**2 + c_hz**2 - 1
    
    # Calculate the discriminant
    discriminant = b**2 - 4 * a * c
    
    if discriminant < 0:
        print("No real solution exists.")
        return
        
    # Calculate the two possible solutions for x
    sqrt_discriminant = math.sqrt(discriminant)
    x1 = (-b + sqrt_discriminant) / (2 * a)
    x2 = (-b - sqrt_discriminant) / (2 * a)
    
    # The geometric relation implies c_hz - c_hb * x >= 0, so x <= c_hz / c_hb.
    # This constraint helps us choose the correct physical solution.
    constraint_val = c_hz / c_hb
    
    solution = None
    if x1 <= constraint_val:
        solution = x1
    if x2 <= constraint_val:
        # If both are valid, there might be an issue. However, given the problem's setup, 
        # only one is expected to be physically meaningful. The other usually arises from squaring.
        # We find that x1 > constraint_val in this specific case.
        solution = x2

    # Let's present the equation with its full numerical coefficients
    c_hb_frac_num = 15
    c_hb_frac_den = 16
    c_hz_frac_num = 9
    c_hz_frac_den = 10

    # Equation: x^2 - (27/16)x + (4409/6400) = 0
    # Or: 6400x^2 - 10800x + 4409 = 0
    print("The problem reduces to solving a quadratic equation for x = lim <b_p, z_p>.")
    print("Assuming h_p, b_p, and z_p are coplanar in the limit, we derive the relation:")
    print("(<h_p, z_p> - <h_p, b_p> * x)^2 = (1 - <h_p, b_p>^2) * (1 - x^2)")
    print(f"Plugging in the values <h_p, b_p> = {c_hb} and <h_p, z_p> = {c_hz}:")
    print(f"({c_hz} - {c_hb}*x)^2 = (1 - {c_hb**2}) * (1 - x^2)")
    print("This simplifies to the quadratic equation:")
    print(f"x^2 - {2 * c_hb * c_hz}*x + {c_hb**2 + c_hz**2 - 1} = 0")
    print("The valid solution for x is:")
    # The expression using fractions is (135 - sqrt(589))/160
    # We output the decimal representation as the final answer.
    print(f"lim <b_p, z_p> = {solution}")
    

solve_for_inner_product()
