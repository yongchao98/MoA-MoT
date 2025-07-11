import math

def count_enclosed_poles():
    """
    Counts the number of poles of the function f(z) inside the contour Gamma.
    """
    pole_count = 0
    # The range for a is given in the problem
    a_min = -2024
    a_max = 2024

    # Determine a reasonable range for k.
    # From analysis, poles must be in the approximate range |u| < 24.
    # |a + 2*pi*k| < 24. Since |a| can be up to 2024,
    # |2*pi*k| < 24 + |a| <= 24 + 2024 = 2048
    # |k| < 2048 / (2*pi) ~= 326
    k_range = 330 

    for a in range(a_min, a_max + 1):
        for k in range(-k_range, k_range + 1):
            # Coordinates of the pole
            x = float(a)
            y = 2 * math.pi * k

            # Transformed coordinates
            u = x + y
            v = x - y

            # The expression for the curve boundary
            # P(u,v) = 1200*v**4 - 20*u**3*v**2 + 3*(u**2 - 400)**3
            
            # Use high precision for calculations
            u2 = u * u
            u3 = u2 * u
            v2 = v * v
            v4 = v2 * v2
            
            term1 = 1200.0 * v4
            term2 = -20.0 * u3 * v2
            
            u2_minus_400 = u2 - 400.0
            term3 = 3.0 * u2_minus_400 * u2_minus_400 * u2_minus_400

            p_val = term1 + term2 + term3

            # Check if the pole is inside the contour
            if p_val < 0:
                pole_count += 1
    
    return pole_count

if __name__ == '__main__':
    N = count_enclosed_poles()
    integral_coeff = 2 * N
    print(f"The equation for the curve Gamma can be written in terms of u=x+y and v=x-y as:")
    print("1200 * v^4 - 20 * u^3 * v^2 + 3 * (u^2 - 400)^3 = 0")
    print(f"\nThe number of poles of f(z) inside Gamma is N = {N}.")
    print(f"The value of the contour integral is 2 * pi * i * N.")
    print(f"So, the integral is 2 * pi * i * {N} = {integral_coeff} * pi * i.")
    # Final answer in the required format
    # The question asks for the value of the integral.
    # The value is an imaginary number. We can't put "i" in the answer format.
    # Usually, when a complex number is requested, one component is requested.
    # Let's assume it's the coefficient of pi*i
    # The format does not allow for letters, so it is likely asking for 2*N
    # print(f"<<<{integral_coeff}>>>")