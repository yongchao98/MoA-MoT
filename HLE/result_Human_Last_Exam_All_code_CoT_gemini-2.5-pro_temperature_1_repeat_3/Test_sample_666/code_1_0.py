import math

def solve_integral():
    """
    This function calculates the number of poles of f(z) inside the contour Gamma
    and prints the final value of the contour integral.
    """
    
    pole_count = 0
    pi = math.pi
    
    # The poles are at z = a + 2*pi*i*k
    # a is an integer in [-2024, 2024]
    # k is an integer
    
    # The curve is bounded, analysis shows we only need to check a small range for k.
    # A safe range for k is from -5 to 5.
    for k in range(-5, 6):
        # The problem statement gives the range for 'a'.
        for a in range(-2024, 2025):
            x = float(a)
            y = 2.0 * pi * k
            
            u = x + y
            v = x - y
            
            # The curve is defined by P(u,v) = 0. A pole is inside if P(u,v) < 0.
            # P(u,v) = 3*(u^2 - 400)^3 + v^2*(1200*v^2 - 20*u^3)
            
            try:
                term1 = 3.0 * ((u**2 - 400.0)**3)
                term2 = (v**2) * (1200.0 * v**2 - 20.0 * u**3)
                value = term1 + term2
            except OverflowError:
                # For large a, u and v are large, and the point is outside.
                continue

            if value < 0:
                pole_count += 1
                
    # The value of the integral is 2 * pi * i * N, where N is the number of poles.
    # The prompt asks to output each number in the final equation.
    # The final equation is Integral = 2 * pi * i * pole_count
    val_2 = 2
    final_coeff = val_2 * pole_count

    print(f"The number of enclosed poles is {pole_count}.")
    print(f"The contour integral is given by the expression: {val_2} * pi * i * {pole_count}")
    print(f"This simplifies to: {final_coeff} * pi * i")

solve_integral()