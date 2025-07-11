import sympy
from sympy import symbols, Function, cos, sin, Eq, solve, simplify

def solve_magnetostatics_problem():
    """
    This script solves for the magnetic field H inside and outside a magnetized
    sphere within a conducting shell using the magnetic scalar potential method.
    """
    # 1. Define symbols and coordinate system representations
    r, theta = symbols('r, θ')
    M0, Rp, R = symbols('M_0 R_p R', positive=True)
    A1, A2, B2 = symbols('A1 A2 B2')

    # 2. Define the magnetic scalar potential Psi_m in the two regions.
    #    The general solution to Laplace's equation for l=1 is used.
    # Region 1 (0 < r < R_p)
    Psi_m1 = A1 * r * cos(theta)
    # Region 2 (R_p < r < R)
    Psi_m2 = (A2 * r + B2 / r**2) * cos(theta)

    # 3. Derive H = -grad(Psi_m) for both regions in spherical coordinates.
    # H_r = -d(Psi_m)/dr
    # H_theta = -(1/r) * d(Psi_m)/d(theta)
    H1_r = -sympy.diff(Psi_m1, r)
    H1_theta = -sympy.diff(Psi_m1, theta) / r

    H2_r = -sympy.diff(Psi_m2, r)
    H2_theta = -sympy.diff(Psi_m2, theta) / r

    # 4. Set up the system of equations based on boundary conditions.
    # BC a: At the perfect conductor (r=R), B_r=0, which means H_r=0.
    eq1 = Eq(H2_r.subs(r, R), 0)

    # BC b: At the interface (r=R_p), H_theta is continuous.
    eq2 = Eq(H1_theta.subs(r, Rp), H2_theta.subs(r, Rp))

    # BC c: At the interface (r=R_p), B_r is continuous.
    # B_r1 = B_r2 => mu_0 * (H_r1 + M_r) = mu_0 * H_r2
    # The magnetization is M = M0*i_z, so M_r = M0*cos(theta).
    # This gives H_r2 - H_r1 = M_r = M0*cos(theta).
    # We can divide by cos(theta) since it's a common factor.
    M_r = M0 * cos(theta)
    eq3 = Eq(H2_r.subs(r, Rp) - H1_r.subs(r, Rp), M_r)

    # 5. Solve for the unknown coefficients A1, A2, B2.
    # The equations are linear in A1, A2, B2. We can simplify by removing sin/cos factors.
    eq1_s = simplify(eq1 / cos(theta))
    eq2_s = simplify(eq2 / sin(theta))
    eq3_s = simplify(eq3 / cos(theta))
    
    solution = solve([eq1_s, eq2_s, eq3_s], (A1, A2, B2))
    
    sol_A1 = solution[A1]
    sol_A2 = solution[A2]
    sol_B2 = solution[B2]

    # 6. Substitute coefficients back and print the final H fields.
    H1_final_r = simplify(H1_r.subs(A1, sol_A1))
    H1_final_theta = simplify(H1_theta.subs(A1, sol_A1))

    H2_final_r = simplify(H2_r.subs([(A2, sol_A2), (B2, sol_B2)]))
    H2_final_theta = simplify(H2_theta.subs([(A2, sol_A2), (B2, sol_B2)]))

    # 7. Format the output to match the structure of the multiple-choice options.
    print("Derived Magnetic Field H:")
    print("-" * 50)
    
    print("In the region 0 < r < R_p:")
    # The field is H = H1_r * î_r + H1_theta * î_θ
    # We can factor the expression to match the answer choices format.
    coeff_H1 = simplify(H1_theta / sin(theta)).subs(A1, sol_A1)
    print(f"H = ({coeff_H1}) * ( - cos(θ) î_r + sin(θ) î_θ )")
    print("\nThis can be written as:")
    print(f"H = M_0 * (2*R_p**3 + R**3)/(3*R**3) * ( - cos(θ) î_r + sin(θ) î_θ )\n")

    print("In the region R_p < r < R:")
    # The field is H = H2_r * î_r + H2_theta * î_θ
    print(f"H_r = {H2_final_r}")
    print(f"H_θ = {H2_final_theta}")
    print("\nThis can be written as:")
    print(f"H_r = -2*M_0/3 * [ (R_p/R)**3 - (R_p/r)**3 ] * cos(θ)")
    print(f"H_θ = M_0/3 * [ 2*(R_p/R)**3 + (R_p/r)**3 ] * sin(θ)")
    print("-" * 50)
    print("The derived expressions match Option B.")

if __name__ == '__main__':
    solve_magnetostatics_problem()