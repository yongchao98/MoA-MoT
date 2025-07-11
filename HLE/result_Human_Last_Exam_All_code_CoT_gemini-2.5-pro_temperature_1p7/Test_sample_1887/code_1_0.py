def solve_set_theory_problem():
    """
    This function formats and prints the solution to the set theory problem.
    The order type of X is omega_(omega+5).
    """
    omega_symbol = "omega"
    subscript_open = "_("
    plus_symbol = "+"
    number_5 = 5
    subscript_close = ")"

    # The final equation is Order Type = omega_(omega+5)
    # The components are printed as requested by the user.
    print("The final order type is: ", end="")
    print(f"{omega_symbol}{subscript_open}{omega_symbol}{plus_symbol}{number_5}{subscript_close}")

solve_set_theory_problem()