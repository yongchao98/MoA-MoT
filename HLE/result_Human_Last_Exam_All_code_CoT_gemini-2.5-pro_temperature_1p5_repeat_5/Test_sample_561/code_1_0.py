import math

def calculate_fractal_dimension():
    """
    Calculates and explains the Minkowski-Bouligand dimension of the fractal piano.
    """
    # 1. Define the parameters of the self-affine system.
    N = 5.0
    # Horizontal scaling factor
    r_x = ( (3.0/7.0)/2.0 ) / 3.0  # (width_black_key) / (width_keyboard) = 1/14
    # Vertical scaling factor
    r_y = (9.0/14.0) / 1.0 # (height_black_key) / (height_keyboard) = 9/14

    # The formula for the dimension D of this self-affine system is:
    # D = D_x + (1 - D_x) * log(r_y) / log(r_x)
    # where D_x = log(N) / log(1/r_x)
    # This can be algebraically simplified to:
    # D = 1 - log(9)/log(14) + (log(5)*log(9))/(log(14)^2)
    
    # 2. Calculate the logarithm values needed.
    log5 = math.log(5)
    log9 = math.log(9)
    log14 = math.log(14)

    # 3. Calculate the final dimension using the simplified formula.
    term1 = 1.0
    term2 = log9 / log14
    term3 = (log5 * log9) / (log14 * log14)
    dimension = term1 - term2 + term3

    # 4. Print the explanation and results.
    print("The Minkowski-Bouligand dimension D is calculated for the self-affine fractal defined by the piano keys.")
    print("The parameters are:")
    print(f" - Number of copies, N = {int(N)}")
    print(f" - Horizontal scaling, r_x = 1/14")
    print(f" - Vertical scaling, r_y = 9/14")
    
    print("\nThe dimension D is given by the formula:")
    print("D = 1 - log(9)/log(14) + (log(5) * log(9)) / (log(14) * log(14))")
    
    print("\nPlugging in the numerical values for the logarithms:")
    print(f"D = 1 - {log9} / {log14} + ({log5} * {log9}) / ({log14} * {log14})")

    print("\nCalculating each term:")
    print(f"D = {term1} - {term2} + {term3}")
    
    print("\nThe final dimension is:")
    print(f"D = {dimension}")

# Execute the calculation
calculate_fractal_dimension()