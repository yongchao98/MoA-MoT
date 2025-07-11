def solve_set_theory_problem():
    """
    Solves the problem by explaining the steps to find the minimal and maximal
    possible cardinalities of X and then calculating their difference.
    """

    print("Step 1: Analyze the given hypotheses.")
    print("The problem is set in ZFC + not CH + 2^omega_1 = omega_3.")
    print("The failure of the continuum hypothesis (not CH) means c = 2^omega_0 > omega_1.")
    print("Since omega_0 < omega_1, we must have c = 2^omega_0 <= 2^omega_1.")
    print("The given hypothesis is 2^omega_1 = omega_3.")
    print("Combining these, we get omega_1 < c <= omega_3.")
    print("This implies that the possible values for the continuum c are omega_2 and omega_3.")
    print("-" * 20)

    print("Step 2: Find the minimal possible cardinality of X.")
    print("To minimize |X|, we seek a model of ZFC where the number of distinct cardinalities of uncountable MAD families is as small as possible.")
    print("It is a known consistency result that Martin's Axiom (MA) + not CH implies that every maximal almost disjoint (MAD) family has cardinality c.")
    print("In such a model, X = {c}, so |X| = 1.")
    print("We need to ensure that MA + not CH is consistent with our hypotheses.")
    print("It is indeed consistent to have a model of ZFC where MA holds, c = omega_2, and 2^omega_1 = omega_3.")
    print("In this model, every uncountable MAD family has size omega_2. Thus, X = {omega_2}.")
    min_cardinality_X = 1
    print(f"The cardinality of X in this model is |X| = {min_cardinality_X}.")
    print("This is the smallest possible non-zero size for X, so the minimal possible cardinality of X is 1.")
    print("-" * 20)

    print("Step 3: Find the maximal possible cardinality of X.")
    print("To maximize |X|, we want a model with the greatest variety of MAD family sizes.")
    print("Powerful consistency results by Saharon Shelah show that it is possible to construct models where MAD families exist for a wide range of cardinalities.")
    print("To maximize the number of possible cardinalities, we should choose the largest possible value for c, which is c = omega_3.")
    print("The uncountable cardinals up to c=omega_3 are omega_1, omega_2, and omega_3.")
    print("All these cardinals (omega_1, omega_2, omega_3) are regular.")
    print("Shelah's results show it is consistent to have a model where a MAD family of size k exists for every regular cardinal k such that omega_1 <= k <= c.")
    print("It is possible to construct such a model which is also consistent with c = omega_3 and 2^omega_1 = omega_3.")
    print("In such a model, there would be MAD families of sizes omega_1, omega_2, and omega_3.")
    print("So, the set X would be {omega_1, omega_2, omega_3}.")
    max_cardinality_X = 3
    print(f"The cardinality of X in this model is |X| = {max_cardinality_X}.")
    print("This represents the maximal possible number of elements for X.")
    print("-" * 20)
    
    print("Step 4: Calculate the difference.")
    difference = max_cardinality_X - min_cardinality_X
    print(f"The maximal possible cardinality of X is {max_cardinality_X}.")
    print(f"The minimal possible cardinality of X is {min_cardinality_X}.")
    print(f"The difference is {max_cardinality_X} - {min_cardinality_X} = {difference}.")

solve_set_theory_problem()