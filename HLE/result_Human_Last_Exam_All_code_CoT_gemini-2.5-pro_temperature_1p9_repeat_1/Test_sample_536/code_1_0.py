import math

def solve_inner_product():
    """
    Calculates the limit of the inner product based on the geometric relationship
    between three coplanar unit vectors.
    """
    # Given limiting inner products
    C_hb = 0.9375
    C_hz = 0.9

    print(f"Given <h, b> = {C_hb}")
    print(f"Given <h, z> = {C_hz}")
    print("\nWe are solving for x = <b, z> using the equation derived from the coplanarity of the vectors:")
    print("x = <h, b> * <h, z> - sqrt((1 - <h, b>^2) * (1 - <h, z>^2))")

    # Calculate the components of the solution
    term1 = C_hb * C_hz
    
    # Calculate squared terms
    C_hb_sq = C_hb**2
    C_hz_sq = C_hz**2
    
    # Calculate terms under the square root
    inside_sqrt1 = 1 - C_hb_sq
    inside_sqrt2 = 1 - C_hz_sq
    
    # The radicand
    radicand = inside_sqrt1 * inside_sqrt2
    
    term2 = math.sqrt(radicand)

    # We choose the smaller root based on the physics of the problem,
    # which corresponds to the minus sign.
    x = term1 - term2

    print("\n--- Calculation Steps ---")
    print(f"First term: <h, b> * <h, z> = {C_hb} * {C_hz} = {term1}")
    print(f"Term under square root: (1 - {C_hb}^2) * (1 - {C_hz}^2) = (1 - {C_hb_sq}) * (1 - {C_hz_sq}) = {radicand}")
    print(f"Second term (sqrt): sqrt({radicand}) = {term2}")
    print("\nFinal Equation:")
    print(f"x = {term1} - {term2}")
    
    print("\n--- Result ---")
    print(f"The value of lim <b_p, z_p> is: {x}")
    
    # Return the value in the specified format
    print(f"\n<<<${x:.4f}$>>>")


solve_inner_product()
