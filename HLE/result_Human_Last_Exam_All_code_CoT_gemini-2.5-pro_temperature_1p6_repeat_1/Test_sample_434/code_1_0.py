def solve_and_print_condition():
    """
    This function prints the derived condition for the initial value y(0)
    that causes the solution of the system of ODEs to blow up, given x(0) > 1.
    """

    # Coefficients and powers in the final equation for the separatrix:
    # y^2 = c1 + c2*x - c3*x^(p_num/p_den)
    # The derived equation is y^2 = 1 + 2*x - 3*x^(2/3)
    c1 = 1
    c2 = 2
    c3 = 3
    p_num = 2
    p_den = 3

    # The solution blows up if y(0) is below the separatrix curve.
    print("For the given system of differential equations and the condition x(0) > 1,")
    print("the solution (x(t), y(t)) blows up in finite time if the initial value y(0) satisfies the following inequality:")
    print("\n---------------------------------------------------------")
    # Using string formatting to output the equation with its numbers
    print(f"y(0) < sqrt( {c1} + {c2}*x(0) - {c3}*x(0)^({p_num}/{p_den}) )")
    print("---------------------------------------------------------\n")
    print("This inequality defines the set of all initial values y(0) (for a given x(0)>1) that result in a blow-up.")

solve_and_print_condition()