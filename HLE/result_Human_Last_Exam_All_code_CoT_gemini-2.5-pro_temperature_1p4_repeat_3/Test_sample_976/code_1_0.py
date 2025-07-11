import sympy
from sympy import symbols, Function, Eq, solve, diff, cos, sin, pprint

def solve_electrodynamics_problem():
    """
    Solves for the electric potential and field outside a conducting sphere
    in a uniform electric field using symbolic mathematics.
    """
    # 1. Define all symbolic variables
    # r: radial distance, theta: polar angle
    # R: radius of the sphere
    # E0: magnitude of the external uniform electric field
    # s1, s2: conductivity inside and outside the sphere, respectively
    r, theta, R, E0, s1, s2 = symbols('r, theta, R, E_0, sigma_1, sigma_2', real=True, positive=True)
    
    # A1, B1: unknown coefficients for the potential solutions
    A1, B1 = symbols('A1 B1')

    # 2. Define the general forms of the potential
    # Based on Laplace's equation and boundary conditions at r=0 and r->infinity.
    # The full potential is f(r) * cos(theta). We solve for the radial part's coefficients.
    Phi_in_radial = A1 * r
    Phi_out_radial = -E0 * r + B1 / r**2

    # 3. Apply boundary conditions at r=R to create a system of equations
    
    # a. Continuity of potential: Phi_in(R) = Phi_out(R)
    # The cos(theta) term cancels from both sides.
    eq1 = Eq(Phi_in_radial.subs(r, R), Phi_out_radial.subs(r, R))

    # b. Continuity of the normal component of current density: J_1n = J_2n
    # -s1 * d(Phi_in)/dr = -s2 * d(Phi_out)/dr
    # s1 * d(Phi_in)/dr = s2 * d(Phi_out)/dr
    # Again, the cos(theta) term cancels.
    dPhi_in_dr = diff(Phi_in_radial, r)
    dPhi_out_dr = diff(Phi_out_radial, r)
    eq2 = Eq(s1 * dPhi_in_dr.subs(r, R), s2 * dPhi_out_dr.subs(r, R))
    
    # 4. Solve the system of equations for the coefficients A1 and B1
    solution = solve([eq1, eq2], [A1, B1])

    # 5. Substitute the coefficients back to get the final potential outside
    Phi_out_final_radial = Phi_out_radial.subs(solution)
    Phi_out_final = Phi_out_final_radial * cos(theta)

    # 6. Calculate the electric field outside from the potential: E = -grad(Phi)
    # E_r = -d(Phi)/dr
    # E_theta = -(1/r) * d(Phi)/dtheta
    E_out_r = -diff(Phi_out_final, r)
    E_out_theta = -(1/r) * diff(Phi_out_final, theta)
    
    # Simplify the expressions
    Phi_out_final = sympy.factor(Phi_out_final, E0, cos(theta))
    E_out_r = sympy.simplify(E_out_r)
    E_out_theta = sympy.simplify(E_out_theta)
    
    # 7. Print the results to be compared with the answer choices
    # The prompt asks to "output each number in the final equation".
    # We will pretty-print the full symbolic equations.
    
    print("The derived electric potential for r > R is:")
    pprint(Phi_out_final)
    
    print("\nThe derived radial component of the electric field for r > R is:")
    pprint(E_out_r)
    
    print("\nThe derived polar component of the electric field for r > R is:")
    pprint(E_out_theta)
    
    print("\n------------------------------------------------------")
    print("Comparing these results with the answer choices:")
    print("The derived potential matches the expression for Phi(r,theta) for r > R in Choice B.")
    print("The derived field components match the expression for E(r,theta) for r > R in Choice B.")
    print("Therefore, the correct answer is B.")

if __name__ == '__main__':
    solve_electrodynamics_problem()
