import cmath

def is_inside(x, y):
    """
    Checks if a point (x, y) is inside the curve Gamma by evaluating the curve's polynomial equation.
    P(x,y) < 0 for points inside the curve.
    """
    # Using float for calculations to handle potential large numbers and non-integer inputs
    x, y = float(x), float(y)
    
    val = (3 * x**6 + 18 * x**5 * y + 45 * x**4 * y**2 + 60 * x**3 * y**3 + 45 * x**2 * y**4 + 18 * x * y**5 + 3 * y**6
           - 20 * x**5 - 20 * x**4 * y + 40 * x**3 * y**2 + 40 * x**2 * y**3 - 20 * x * y**4 - 20 * y**5
           - 2400 * x**4 - 19200 * x**3 * y - 14400 * x**2 * y**2 - 19200 * x * y**3 - 2400 * y**4
           + 1440000 * x**2 + 2880000 * x * y + 1440000 * y**2
           - 192000000)
    
    return val < 0

def solve_contour_integral():
    """
    Calculates the contour integral by counting the poles of f(z) inside Gamma.
    """
    pole_count = 0
    
    # The curve is bounded, so we only need to check a finite range of 'a' and 'k'.
    # Analysis shows |y| < 15, so |2*pi*k| < 15 => k is in {-2, -1, 0, 1, 2}.
    # Analysis shows |x| < 20, so we check 'a' in a slightly larger range for safety.
    a_range = range(-25, 26) 
    k_range = range(-3, 4)

    for a in a_range:
        # The sum for f(z) is over a from -2024 to 2024.
        if not (-2024 <= a <= 2024):
            continue
        for k in k_range:
            pole_x = a
            pole_y = 2 * cmath.pi * k
            
            # cmath.pi is complex, take the real part of y for the check
            if is_inside(pole_x, pole_y.real):
                pole_count += 1
                
    # The value of the integral is 2 * pi * i * N, where N is the number of poles.
    integral_value = 2 * cmath.pi * 1j * pole_count
    
    print(f"The number of poles inside the curve is N = {pole_count}")
    print("The value of the contour integral is given by the formula: 2 * pi * i * N")
    print(f"Using the values: 2 * {cmath.pi} * {1j} * {pole_count}")
    print(f"The final result is: {integral_value}")

solve_contour_integral()