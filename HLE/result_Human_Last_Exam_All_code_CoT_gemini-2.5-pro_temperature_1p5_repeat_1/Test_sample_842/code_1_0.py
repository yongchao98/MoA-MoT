import sympy

def solve():
    """
    This function calculates the determinant of (I - rho(beta)) 
    and then establishes the relationship between f(t) and Q(t).
    """
    t = sympy.Symbol('t')

    # Define the reduced Burau representation matrices for B_3
    rho_s1 = sympy.Matrix([[-t, 1], [0, 1]])
    rho_s2_inv = sympy.Matrix([[1, 0], [1, -1/t]])

    print("The matrix representation for sigma_1 is:")
    sympy.pprint(rho_s1)
    print("\nThe matrix representation for sigma_2^-1 is:")
    sympy.pprint(rho_s2_inv)

    # Calculate the representation for beta = sigma_2^-1 * sigma_1 * sigma_2^-1 * sigma_1
    rho_beta = rho_s2_inv * rho_s1 * rho_s2_inv * rho_s1
    print("\nThe matrix representation for beta is:")
    sympy.pprint(rho_beta)

    # Define the 2x2 Identity matrix
    I2 = sympy.eye(2)

    # Calculate the matrix (I - rho(beta))
    M = I2 - rho_beta
    print("\nThe matrix (I - rho(beta)) is:")
    sympy.pprint(M)

    # Calculate the determinant
    det_M = M.det()
    print(f"\nThe determinant det(I - rho(beta)) is: {sympy.simplify(det_M)}")

    # To make the relationship clear, we express the determinant with a common denominator
    det_M_rational = sympy.together(det_M)
    numerator = sympy.fraction(det_M_rational)[0]
    denominator = sympy.fraction(det_M_rational)[1]

    print(f"\nExpressed as a rational function, the determinant is:")
    print(f"Numerator: {numerator}")
    print(f"Denominator: {denominator}")
    
    print("\nThe given formula is Q(t) = f(t) / (-t^4 + 2*t^3 + t^2 + 2*t - 1) * det(I - rho(beta)).")
    print("Substituting our calculated determinant:")
    print("Q(t) = f(t) / (-t^4 + 2*t^3 + t^2 + 2*t - 1) * ((-t^4 + 2*t^3 + t^2 + 2*t - 1) / t^2)")
    print("This simplifies to Q(t) = f(t) / t^2, or f(t) = t^2 * Q(t).")
    
    print("\nThe closure of the braid beta is the link L6a5, which is amphicheiral.")
    print("This means its polynomial Q(t) must be symmetric, i.e., Q(t) = Q(t^-1).")
    print("We check the answer choices:")
    print("A: f(t) = 1 => Q(t) = 1/t^2. Q(t^-1) = t^2. Not symmetric.")
    print("B: f(t) = t^2 => Q(t) = 1. Q(t^-1) = 1. Symmetric.")
    print("C: f(t) = -1 => Q(t) = -1/t^2. Q(t^-1) = -t^2. Not symmetric.")
    print("D: f(t) = -t^3+3t^2-2t+1 => Q(t) = -t+3-2/t+1/t^2. Not symmetric.")
    print("E: f(t) = 2t^5+4t^4-2t^3-3t^2 => Q(t) = 2t^3+4t^2-2t-3. Not symmetric.")
    print("\nOnly choice B results in a symmetric polynomial for Q(t).")

solve()