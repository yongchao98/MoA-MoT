import sympy as sp

def solve_and_verify():
    """
    This function defines the potential expressions from Option A and verifies
    that they satisfy all boundary conditions of the problem using symbolic math.
    """
    # Define symbols
    x, y, k, a, b, sigma0, eps1, eps2 = sp.symbols(
        'x y k a b sigma_0 epsilon_1 epsilon_2', real=True, positive=True
    )
    sp.init_printing(use_unicode=False, wrap_line=False)

    # Surface charge density as given in the problem
    sigma_s = sigma0 * sp.sin(k*x)

    # We test Option A, which was found through manual derivation.
    # Denominator from Option A
    denominator = k * (eps2 * sp.cosh(k*a) * sp.sinh(k*b) + eps1 * sp.sinh(k*a) * sp.cosh(k*b))

    # Potential in Region 2 (0 < y < a) from Option A
    phi2_numerator = -sigma0 * sp.sinh(k*b) * sp.sinh(k*(y - a)) * sp.sin(k*x)
    Phi2 = phi2_numerator / denominator

    # Potential in Region 1 (-b < y < 0) from Option A
    phi1_numerator = sigma0 * sp.sinh(k*a) * sp.sinh(k*(y + b)) * sp.sin(k*x)
    Phi1 = phi1_numerator / denominator

    print("--- Verifying Option A against the Boundary Conditions ---")

    # BC1: Potential at y = -b must be 0
    check1 = sp.simplify(Phi1.subs(y, -b))
    print(f"1. Potential at y = -b: Phi_1(y=-b) = {check1} (Success if 0)")

    # BC2: Potential at y = a must be 0
    check2 = sp.simplify(Phi2.subs(y, a))
    print(f"2. Potential at y = a: Phi_2(y=a) = {check2} (Success if 0)")

    # BC3: Potential must be continuous at y = 0
    phi1_at_0 = Phi1.subs(y, 0)
    phi2_at_0 = Phi2.subs(y, 0)
    check3 = sp.simplify(phi1_at_0 - phi2_at_0)
    print(f"3. Continuity at y = 0: simplify(Phi_1(0) - Phi_2(0)) = {check3} (Success if 0)")

    # BC4: Gauss's Law at the interface y=0
    # The condition is: eps2 * d(Phi2)/dy - eps1 * d(Phi1)/dy = -sigma_s
    dPhi1_dy = sp.diff(Phi1, y)
    dPhi2_dy = sp.diff(Phi2, y)
    
    dPhi1_dy_at_0 = dPhi1_dy.subs(y, 0)
    dPhi2_dy_at_0 = dPhi2_dy.subs(y, 0)
    
    # Check if the jump condition holds by seeing if the expression simplifies to 0
    jump_check = sp.simplify(eps2 * dPhi2_dy_at_0 - eps1 * dPhi1_dy_at_0 + sigma_s)
    print(f"4. Gauss's Law at y=0: simplify(eps2*d(Phi2)/dy - eps1*d(Phi1)/dy + sigma_s) = {jump_check} (Success if 0)")
    
    if check1 == 0 and check2 == 0 and check3 == 0 and jump_check == 0:
        print("\nConclusion: Option A satisfies all boundary conditions.")
    else:
        print("\nConclusion: Option A does NOT satisfy all boundary conditions.")

    # The question asks for the potential in the region 0 <= y <= a.
    # We print this expression, including all symbolic variables as requested.
    print("\nThe electric potential Phi(x, y) in the region 0 <= y <= a is:")
    phi2_string_num = f"-sigma_0 * sinh(k*b) * sinh(k*(y - a)) * sin(k*x)"
    phi2_string_den = f"k * (epsilon_2*cosh(k*a)*sinh(k*b) + epsilon_1*sinh(k*a)*cosh(k*b))"
    
    # We print each "number" (symbol) in the equation by constructing the formula string.
    print(f"    Phi(x, y) = ({phi2_string_num}) / ({phi2_string_den})")

solve_and_verify()