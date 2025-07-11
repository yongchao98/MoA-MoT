def solve_scl():
    """
    Computes the stable commutator length of the element c in the group G.

    Let F_i be the free group with basis {a_i, b_i}. Let c_i = [a_i, b_i].
    Let G be the free product of F_i for i = 1, ..., 19.
    Let c be the product c = c_1^30 * c_2^30 * ... * c_19^30.

    We compute the stable commutator length (scl) of c in G.
    """

    # Number of free groups in the free product
    num_groups = 19

    # The exponent of each commutator
    exponent = 30

    print("Problem setup:")
    print(f"Group G = F_1 * F_2 * ... * F_{num_groups}")
    print(f"Element c = c_1^{exponent} * c_2^{exponent} * ... * c_{num_groups}^{exponent}")
    print("-" * 30)

    # Step 1: Apply Fujiwara's formula for scl in free products.
    # scl_G(g) = (m - sum(||g_k||_scl)) / 2
    # where g = g_1 * ... * g_m is a cyclically reduced word.
    print("Applying Fujiwara's formula for stable commutator length (scl):")
    print("scl_G(c) = (m - sum_of_norms) / 2\n")

    # Step 2: Identify the parameters in the formula.
    # The element c is a product of num_groups syllables, g_i = c_i^exponent.
    # This word is cyclically reduced.
    m = num_groups
    print(f"The number of syllables in the word c is m = {m}.")

    # Step 3: Calculate the scl norm of each syllable.
    # The norm ||g_i||_scl is the scl norm of the homology class of g_i in H_1(F_i, R).
    # Each syllable g_i = c_i^exponent = [a_i, b_i]^exponent is in the commutator subgroup [F_i, F_i].
    # Therefore, its image in the homology group H_1(F_i, R) is the zero element.
    # The scl norm of the zero element is 0.
    norm_of_each_syllable = 0
    print(f"Each syllable c_i^{exponent} is a commutator, so its homology class is 0.")
    print(f"The scl norm of the 0 homology class is 0.")

    # The sum of the norms is the sum of num_groups zeros.
    sum_of_norms = num_groups * norm_of_each_syllable
    print(f"The sum of the norms of all syllables is {sum_of_norms}.")
    print("-" * 30)

    # Step 4: Compute the final scl value.
    scl_c = (m - sum_of_norms) / 2.0

    print("Final Calculation:")
    # The user requested to see the numbers in the final equation.
    print(f"scl(c) = ({m} - {sum_of_norms}) / 2")
    print(f"scl(c) = {scl_c}")

solve_scl()