def solve_auditor_choice():
    """
    This function determines and prints the set of values for theta for which
    a company would choose a lenient auditor over a strict one.

    The logic is as follows:
    1.  Establish the cost functions for choosing a lenient (x=0) and a strict (x=1) auditor.
        - C(theta, 0): Cost with a lenient auditor.
        - C(theta, 1): Cost with a strict auditor.
    2.  A firm chooses the lenient auditor if C(theta, 0) < C(theta, 1).
    3.  Derive the inequality that represents this choice.
    4.  Print the resulting inequality.
    """

    # The problem asks to output the numbers in the final equation.
    # The final inequality is t(theta) < p(1) * F(theta).
    # The number in this equation is 1, which is the index for the strict auditor.
    strict_auditor_index = 1

    print("Based on the economic model, a company of type theta chooses a lenient auditor if its tax liability is less than the expected penalty from hiring a strict auditor.")
    print("The set of values for theta for which a company chooses a lenient auditor is defined by the following inequality:")
    
    # We use an f-string to embed the number into the printed equation.
    final_equation = f"t(theta) < p({strict_auditor_index}) * F(theta)"
    print(final_equation)

solve_auditor_choice()