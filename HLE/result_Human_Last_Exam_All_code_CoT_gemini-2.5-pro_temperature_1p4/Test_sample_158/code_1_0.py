import sympy
from sympy import symbols, Eq, solve, pretty_print, factor

def solve_magnetic_shell_problem():
    """
    Solves for the permeability of a cylindrical shell and the interior magnetic
    field under the condition that an external uniform magnetic field is not
    distorted outside the shell. The derivation is performed symbolically.
    """
    print("This script solves a magnetostatics problem for a cylindrical shell.")
    print("The goal is to find the shell's permeability and the interior magnetic field.")
    print("-" * 70)

    # 1. Define symbolic variables
    H0, R1, R2, mu0 = symbols('H_0 R_1 R_2 mu_0', real=True, positive=True)
    mu_r, A1, A2, B2 = symbols('mu_r A_1 A_2 B_2', real=True)

    print("Step 1: Define the magnetic scalar potential Phi_m in each region.")
    print("The applied field H_ext = H_0 * x_hat corresponds to Phi_m_ext = -H_0 * rho * cos(phi).")
    print("\nPotentials (proportional to cos(phi)):")
    print(f"  Interior (rho < R1):      Phi_1 = A_1 * rho")
    print(f"  Shell (R1 < rho < R2):    Phi_2 = A_2 * rho + B_2/rho")
    print(f"  Exterior (rho > R2):      Phi_3 = -H_0 * rho  (undistorted field)")
    print("(The common cos(phi) factor is omitted for simplicity in the equations)")
    print("-" * 70)

    # 2. Set up boundary condition equations
    print("Step 2: Set up boundary condition equations at rho=R1 and rho=R2.")
    print("BC 1: Tangential H is continuous -> Phi_m is continuous.")
    print("BC 2: Normal B is continuous -> mu * d(Phi_m)/d(rho) is continuous.")
    
    # Boundary at rho = R2
    # Eq1: Phi_2(R2) = Phi_3(R2)
    eq1 = Eq(A2 * R2 + B2 / R2, -H0 * R2)
    # Eq2: mu * d(Phi_2)/d(rho)|R2 = mu0 * d(Phi_3)/d(rho)|R2  (where mu = mu_r*mu0)
    # => mu_r * (A2 - B2/R2**2) = 1 * (-H0)
    eq2 = Eq(mu_r * (A2 - B2 / R2**2), -H0)
    
    # Boundary at rho = R1
    # Eq3: Phi_1(R1) = Phi_2(R1)
    eq3 = Eq(A1 * R1, A2 * R1 + B2 / R1)
    # Eq4: mu0 * d(Phi_1)/d(rho)|R1 = mu * d(Phi_2)/d(rho)|R1
    # => 1 * A1 = mu_r * (A2 - B2/R1**2)
    eq4 = Eq(A1, mu_r * (A2 - B2 / R1**2))
    
    print("\nSystem of equations from boundary conditions:")
    pretty_print(eq1)
    pretty_print(eq2)
    pretty_print(eq3)
    pretty_print(eq4)
    print("-" * 70)

    # 3. Solve for coefficients A2 and B2
    print("Step 3: Solve for A2 and B2 in terms of mu_r from the R2 boundary equations.")
    sol_A2_B2 = solve([eq1, eq2], (A2, B2))
    A2_sol = sol_A2_B2[A2]
    B2_sol = sol_A2_B2[B2]
    print("Solution for A2:")
    pretty_print(Eq(A2, A2_sol))
    print("\nSolution for B2:")
    pretty_print(Eq(B2, B2_sol))
    print("-" * 70)

    # 4. Find the condition on mu_r
    print("Step 4: Establish a consistency condition from the R1 boundary equations.")
    # Substitute the expressions for A1 from eq3 and eq4 into each other.
    # A2*R1 + B2/R1 = R1 * mu_r * (A2 - B2/R1**2)
    # A2 + B2/R1**2 = mu_r * (A2 - B2/R1**2)
    consistency_eq = Eq(A2 + B2 / R1**2, mu_r * (A2 - B2 / R1**2))
    print("This condition is: ")
    pretty_print(consistency_eq)
    
    print("\nSubstituting A2 and B2 and simplifying leads to a condition on mu_r:")
    # A2(1-mu_r) + B2/R1**2 * (1+mu_r) = 0
    final_cond_eq = A2_sol * (1 - mu_r) + (B2_sol / R1**2) * (1 + mu_r)
    # The simplification gives (mu_r**2 - 1) multiplied by a non-zero term.
    # So we solve for mu_r**2 - 1 = 0
    mu_r_solutions = solve(factor(final_cond_eq), mu_r)
    mu_r_nontrivial = mu_r_solutions[1]

    print("\n======================= FINAL RESULTS =======================")
    print("\n1. Required Permeability of the Shell")
    print("Solving the consistency equation for the relative permeability gives two solutions:")
    print(f"  mu_r = {mu_r_solutions[0]} (trivial case, shell is a vacuum)")
    print(f"  mu_r = {mu_r_solutions[1]} (the required non-trivial solution)")

    final_mu_expr = mu_r_nontrivial * mu0
    print("\nThe required permeability of the shell is mu = -1 * mu_0.")
    print("This value corresponds to an idealized 'perfect diamagnet'.")
    print("\nFinal Equation for Permeability:")
    pretty_print(Eq(symbols('mu'), final_mu_expr))
    print("The number in this equation is -1.")

    print("\n2. Magnetic Field in the Interior Region")
    # Substitute mu_r = -1 back into the solutions for A2 and B2
    # The limit must be taken as sympy will have 0 in the denominator
    A2_final = sympy.limit(A2_sol, mu_r, mu_r_nontrivial)
    B2_final = sympy.limit(B2_sol, mu_r, mu_r_nontrivial)

    # Now find A1 from Eq3
    A1_final = sympy.simplify((A2_final * R1 + B2_final / R1) / R1)
    
    print(f"Using mu_r = {mu_r_nontrivial}, the coefficient A1 for the interior potential is found to be:")
    pretty_print(Eq(A1, A1_final))
    
    print("\nThe interior magnetic field H_int is uniform and given by H_int = -A1 * x_hat.")
    print("\nFinal Equation for the Interior Magnetic Field:")
    final_H_int_expr = -A1_final * symbols('x_hat')
    pretty_print(Eq(symbols('H_int'), final_H_int_expr))
    print("The numbers and symbols in the final equation for the magnitude of H_int are:")
    print("the external field magnitude H_0, the outer radius R_2 raised to the power of 2,")
    print("and the inner radius R_1 raised to the power of 2.")

if __name__ == '__main__':
    solve_magnetic_shell_problem()