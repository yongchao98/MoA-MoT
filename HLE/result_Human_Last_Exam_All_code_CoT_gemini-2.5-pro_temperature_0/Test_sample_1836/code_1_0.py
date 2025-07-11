def solve_ordinal_problem():
    """
    This script explains the step-by-step solution to the given set theory problem.
    It determines the order type of the set Y and answers the final question.
    """

    print("--- Step 1: Characterizing the sets kappa_n ---")
    print("Let kappa be a measurable cardinal. We define a sequence of sets of ordinals.")
    print("kappa_0 = kappa = {alpha | alpha < kappa}")
    print("kappa_1 is the set of successor ordinals in kappa_0. A successor ordinal is of the form beta + 1.")
    print("So, kappa_1 = {beta + 1 | beta < kappa}.")
    print("kappa_2 is the set of successor ordinals in kappa_1. The successor of an element (beta + 1) in kappa_1 is (beta + 1) + 1 = beta + 2.")
    print("So, kappa_2 = {beta + 2 | beta < kappa}.")
    print("By induction, we can establish the general form for any n in omega = {0, 1, 2, ...}:")
    print("kappa_n = {beta + n | beta < kappa}\n")

    print("--- Step 2: Determining the intersection Y ---")
    print("Y is the intersection of all kappa_n for n in omega.")
    print("An ordinal alpha is in Y if and only if alpha is in kappa_n for all n.")
    print("1. If alpha is in Y, it must be in kappa_0, so alpha < kappa.")
    print("2. For any n >= 1, alpha must be in kappa_n. This means alpha can be written as beta_n + n for some beta_n < kappa.")
    print("3. The condition alpha = beta_n + n implies that alpha >= n for all n in omega.")
    print("4. Therefore, alpha must be greater than or equal to the supremum of omega, which is omega itself (alpha >= omega).")
    print("Combining these, an ordinal alpha is in Y if and only if: omega <= alpha < kappa.")
    print("So, Y = {alpha | omega <= alpha < kappa}\n")

    print("--- Step 3: Calculating the order type of Y ---")
    print("The order type of Y, ot(Y), is the 'length' of the set of ordinals from omega up to (but not including) kappa.")
    print("In ordinal arithmetic, this is written as: ot(Y) = kappa - omega.")
    print("Ordinal subtraction (beta - alpha) finds the unique ordinal gamma such that alpha + gamma = beta.")
    print("We need to find gamma such that: omega + gamma = kappa.")
    print("For any infinite cardinal kappa (which a measurable cardinal is), and any ordinal alpha < kappa, it is a standard result that alpha + kappa = kappa.")
    print("Setting gamma = kappa, we get omega + kappa = kappa. This is the correct unique solution.")
    print("Therefore, ot(Y) = kappa.\n")

    print("--- Step 4: Answering the final question ---")
    print("The question is: For how many ordinals alpha is the order type of Y at least alpha?")
    print("We found ot(Y) = kappa. So the question is: For how many ordinals alpha is kappa >= alpha?")
    print("The set of ordinals satisfying this is S = {alpha | alpha <= kappa}.")
    print("This set is well-ordered. In set theory, the most complete answer for 'how many' ordinals satisfy a property is the order type of the set of those ordinals.")
    print("The set S contains all ordinals up to kappa, which has an order type of kappa + 1.\n")

    print("--- Final Answer ---")
    print("The number of ordinals alpha for which ot(Y) >= alpha is given by the order type of the set {alpha | alpha <= kappa}.")
    print("The final equation for this order type is:")
    
    part1 = "kappa"
    part2 = "1"
    print(f"Order Type = {part1} + {part2}")

solve_ordinal_problem()