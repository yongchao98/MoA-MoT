import math

def solve_contour_integral():
    """
    This function calculates the specified contour integral by:
    1. Defining the simplified equation of the curve Gamma.
    2. Counting the number of poles of the function f(z) that lie inside Gamma.
    3. Calculating the integral using the Residue Theorem.
    """

    # The simplified equation of the curve is G(p,q) = 0, where p = x+y and q = x-y.
    # The point (x,y) is inside the curve if G(p,q) < 0.
    def G(p, q):
        # The coefficients of the simplified polynomial G(p,q) are:
        c_p6 = 3
        c_p3q2 = -20
        c_p4 = -3600
        c_q4 = 1200
        c_p2 = 1440000
        c_const = -192000000
        
        # This is G(p,q):
        # 3*p^6 - 20*p^3*q^2 - 3600*p^4 + 1200*q^4 + 1440000*p^2 - 192000000
        val = (c_p6 * p**6 + c_p3q2 * p**3 * q**2 + c_p4 * p**4 + 
               c_q4 * q**4 + c_p2 * p**2 + c_const)
        return val

    # Based on an analysis of the curve's equation, the search space for poles
    # z = a + 2*pi*i*n can be bounded.
    # The range for 'a' is integers from -17 to 17.
    # The range for 'n' is integers from -2 to 2.
    a_range = range(-17, 18)
    n_range = range(-2, 3)

    pole_count = 0
    for a in a_range:
        for n in n_range:
            # The coordinates of the pole z = x + iy
            x = float(a)
            y = 2 * math.pi * n
            
            # The transformed coordinates
            p = x + y
            q = x - y
            
            # Check if the pole (x,y) is inside the curve Gamma
            if G(p, q) < 0:
                pole_count += 1

    # By the Residue Theorem, the integral is 2 * pi * i * (sum of residues).
    # Each pole has a residue of 1.
    # So the integral is 2 * pi * i * N, where N is the number of poles.
    
    # "Final equation" components
    factor_2 = 2
    factor_pi = math.pi
    factor_i = 1j
    num_poles = pole_count

    integral_value = factor_2 * factor_pi * factor_i * num_poles
    
    # Print the components of the calculation as requested
    print(f"The number of poles inside the curve Gamma is N = {num_poles}.")
    print("The value of the contour integral is given by the equation: 2 * pi * i * N.")
    print(f"Plugging in the numbers: {factor_2} * {factor_pi} * i * {num_poles} = {integral_value}")

solve_contour_integral()