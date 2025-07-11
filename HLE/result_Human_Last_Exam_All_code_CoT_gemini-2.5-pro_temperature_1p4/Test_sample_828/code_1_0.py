import sympy as sp

def solve_for_K():
    """
    This function performs the coordinate transformation of the stress-energy tensor
    and solves for the unknown K.
    """
    # Step 1: Define all the necessary mathematical symbols
    # T_cal represents the calligraphic T (mathcal{T})
    T_cal, a, omega, theta, r = sp.symbols('mathcal{T} a omega theta r')
    K = sp.Symbol('K')

    # Step 2: Define the non-zero components of the tensor in Cartesian coordinates
    # These are given as functions of the polar angle theta
    T_xx = T_cal * (a * omega)**2 * sp.sin(theta)**2
    T_yy = T_cal * (a * omega)**2 * sp.cos(theta)**2
    T_xy = -T_cal * (a * omega)**2 * sp.sin(theta) * sp.cos(theta)
    
    # Step 3: Define the partial derivatives needed for the transformation
    # The transformation is x = r*cos(theta), y = r*sin(theta)
    dx_dtheta = sp.diff(r * sp.cos(theta), theta)
    dy_dtheta = sp.diff(r * sp.sin(theta), theta)

    # Step 4: Apply the tensor transformation law to find the T_thetatheta component in polar coordinates
    # T'_thetatheta = (dx/dtheta)^2*T_xx + (dy/dtheta)^2*T_yy + 2*(dx/dtheta)*(dy/dtheta)*T_xy
    T_thetatheta_prime = (dx_dtheta)**2 * T_xx + \
                         (dy_dtheta)**2 * T_yy + \
                         2 * dx_dtheta * dy_dtheta * T_xy

    # The problem describes a ring of radius 'a', so we evaluate the transformed
    # component at r = a
    T_thetatheta_prime_on_ring = T_thetatheta_prime.subs(r, a)
    
    # Simplify the resulting expression
    T_thetatheta_prime_simplified = sp.simplify(T_thetatheta_prime_on_ring)

    # Step 5: The problem posits that T_thetatheta = a^2 * sin^2(theta) * T + K
    # We solve for K by equating our result with the problem's statement
    # K = T_thetatheta' - (a^2 * sin^2(theta) * T)
    K_solution = T_thetatheta_prime_simplified - (a**2 * sp.sin(theta)**2 * T_cal)
    
    # Step 6: Print the components of the final equation for K as requested.
    # The equation for K is: K = a**4 * omega**2 * T_cal - a**2 * sin(theta)**2 * T_cal
    
    # Extract coefficients and terms to print them individually
    # K = (1)*(a**4 * omega**2 * T_cal) + (-1)*(a**2 * sin(theta)**2 * T_cal)
    term1 = T_thetatheta_prime_simplified
    term2 = -(a**2 * sp.sin(theta)**2 * T_cal)
    
    print("The final equation for K is composed of two terms:")
    print(f"Term 1: {term1}")
    print(f"Term 2: {term2}")
    
    # Re-writing with explicit numbers as requested by the prompt
    print("\nIn the format 'number * expression', the final equation is:")
    print(f"K = (1) * {term1} + (-1) * {sp.Abs(term2)}")

    print(f"\nThus, the expression for K is:")
    print(f"K = {K_solution}")

if __name__ == '__main__':
    solve_for_K()