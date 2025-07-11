import sympy

def solve_steady_state():
    """
    This function derives and prints the formula for the steady-state probability pi_0
    for the given birth-death process.
    """
    # Define rho as a symbolic variable
    rho = sympy.Symbol('rho')

    # The sum of the ratios pi_n / pi_0 from n=0 to infinity is the sum of rho^n / n!
    # which is the Taylor series for e**rho.
    # sum(pi_n) = pi_0 * sum(rho**n / n!) = 1
    # pi_0 * e**rho = 1
    # pi_0 = 1 / e**rho = e**(-rho)
    
    # Define pi_0 as a function of rho
    pi_0 = sympy.exp(-rho)
    
    # Create the final equation as a symbolic expression
    equation = sympy.Eq(sympy.Symbol('pi_0'), pi_0)
    
    # Print the final result
    # The format "output each number in the final equation" is interpreted as
    # printing the components of the derived symbolic formula.
    print("The derived steady-state probability pi_0 is:")
    print(f"{equation.lhs} = {equation.rhs}")

if __name__ == "__main__":
    solve_steady_state()