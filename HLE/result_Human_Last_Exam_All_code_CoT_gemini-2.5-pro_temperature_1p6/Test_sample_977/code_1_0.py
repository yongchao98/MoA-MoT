import sympy
from sympy import sin, cos, sinh, cosh, diff, simplify, Eq

def solve_potential():
    """
    Solves the electrostatics problem by verifying the answer choices against the boundary conditions.
    """
    # Define the symbolic variables
    x, y, k, a, b = sympy.symbols('x y k a b', real=True, positive=True)
    sigma_0, eps1, eps2 = sympy.symbols('sigma_0 epsilon_1 epsilon_2', real=True, positive=True)

    # Define the surface charge distribution
    sigma_s = sigma_0 * sin(k*x)

    # Dictionary to store the potential functions for each answer choice
    # Each entry contains a tuple (Phi_1 for -b<y<0, Phi_2 for 0<y<a)
    potentials = {}

    # Common denominator for choices A, D, E
    Den_A = k * (eps2 * cosh(k*a) * sinh(k*b) + eps1 * sinh(k*a) * cosh(k*b))
    
    # Choice A
    Phi2_A = (-sigma_0 * sinh(k*b) * sinh(k*(y - a)) * sin(k*x)) / Den_A
    Phi1_A = (sigma_0 * sinh(k*a) * sinh(k*(y + b)) * sin(k*x)) / Den_A
    potentials['A'] = (Phi1_A, Phi2_A)

    # Choice B
    Den_B = k * (eps2 * sinh(k*a) * cosh(k*b) + eps1 * cosh(k*a) * sinh(k*b))
    Phi2_B = (-sigma_0 * sinh(k*a) * sinh(k*(y - a)) * sin(k*x)) / Den_B
    Phi1_B = (sigma_0 * sinh(k*b) * sinh(k*(y + b)) * sin(k*x)) / Den_B
    potentials['B'] = (Phi1_B, Phi2_B)

    # Choice C
    Den_C = k * (eps2 * cosh(k*a) * cosh(k*b) + eps1 * sinh(k*a) * sinh(k*b))
    Phi2_C = (-sigma_0 * sinh(k*b) * sinh(k*(y - a)) * sin(k*x)) / Den_C
    Phi1_C = (sigma_0 * sinh(k*a) * sinh(k*(y + b)) * sin(k*x)) / Den_C
    potentials['C'] = (Phi1_C, Phi2_C)
    
    # Choice D
    Phi2_D = (sigma_0 * sinh(k*b) * sinh(k*(y - a)) * sin(k*x)) / Den_A
    Phi1_D = (-sigma_0 * sinh(k*a) * sinh(k*(y + b)) * sin(k*x)) / Den_A
    potentials['D'] = (Phi1_D, Phi2_D)
    
    # Choice E
    Phi2_E = (sigma_0 * sinh(k*b) * sinh(k*(y-a)) * sin(k*x)) / Den_A
    Phi1_E = (sigma_0 * sinh(k*a) * sinh(k*(y+b)) * sin(k*x)) / Den_A
    potentials['E'] = (Phi1_E, Phi1_A)


    # Check each answer choice
    correct_answer = None
    for choice, (phi1, phi2) in potentials.items():
        print(f"--- Checking Choice {choice} ---")

        # BC1: Potential is zero at y = a
        bc1_check = simplify(phi2.subs(y, a)) == 0
        print(f"BC1 (Phi(y=a)=0) satisfied: {bc1_check}")

        # BC2: Potential is zero at y = -b
        bc2_check = simplify(phi1.subs(y, -b)) == 0
        print(f"BC2 (Phi(y=-b)=0) satisfied: {bc2_check}")

        # BC3: Potential is continuous at y = 0
        phi1_at_0 = phi1.subs(y, 0)
        phi2_at_0 = phi2.subs(y, 0)
        bc3_check = simplify(phi1_at_0 - phi2_at_0) == 0
        print(f"BC3 (Continuity at y=0) satisfied: {bc3_check}")

        # BC4: Discontinuity in Electric Displacement Field at y = 0
        d_phi1_dy = diff(phi1, y).subs(y, 0)
        d_phi2_dy = diff(phi2, y).subs(y, 0)
        displacement_discontinuity = eps1 * d_phi1_dy - eps2 * d_phi2_dy
        bc4_check = simplify(displacement_discontinuity - sigma_s) == 0
        print(f"BC4 (E-field discontinuity at y=0) satisfied: {bc4_check}")

        if all([bc1_check, bc2_check, bc3_check, bc4_check]):
            correct_answer = choice
            print(f"\nChoice {choice} is the correct solution.")
            print("\nThe electric potential Phi(x, y) in the region 0 <= y <= a is:")
            
            # For better readability, substitute long names
            s_0_sym = sympy.Symbol('sigma_0')
            k_sym = sympy.Symbol('k')
            a_sym = sympy.Symbol('a')
            b_sym = sympy.Symbol('b')
            eps_1_sym = sympy.Symbol('epsilon_1')
            eps_2_sym = sympy.Symbol('epsilon_2')
            
            pretty_phi2 = phi2.subs({sigma_0: s_0_sym, k:k_sym, a:a_sym, b:b_sym, eps1:eps_1_sym, eps2:eps_2_sym})
            
            print(f"Phi(x, y) = {sympy.pretty(pretty_phi2, use_unicode=False)}")

    return correct_answer

if __name__ == '__main__':
    final_answer = solve_potential()
    if final_answer:
        print(f"\n<<<A>>>")
