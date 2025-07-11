def solve_sqs_doubling_question():
    """
    This script addresses a three-part question about the SQS doubling construction
    by explaining the reasoning and providing calculations for the specific case of v=4.
    """

    # Part (a): True or False: Each element is contained in exactly v-1 ND-pairs.
    v_a = 4
    # The number of blocks (and thus primary ND-pairs) containing any point in an SQS(2v)
    # is r' = (2v - 1)(2v - 2) / 6.
    r_prime = (2 * v_a - 1) * (2 * v_a - 2) / 6
    v_minus_1_a = v_a - 1
    
    print("(a) False")
    print("The number of ND-pairs containing any single point in a nested SQS(2v) is given by the formula r' = (2v-1)(2v-2)/6.")
    print(f"For v = {v_a}, r' = (2*{v_a}-1)*(2*{v_a}-2)/6 = {int(r_prime)}.")
    print(f"The question asks if this equals v-1, which is {v_minus_1_a} for v = {v_a}.")
    print(f"Since {int(r_prime)} != {v_minus_1_a}, the statement is false for v >= 4.")
    print("-" * 20)

    # Part (b): Multiplicity relationship for horizontal pairs.
    print("(b) \u03BC") # \u03BC is unicode for mu
    print("A standard doubling construction relates the multiplicity of a horizontal pair {(x, 0), (y, 0)} directly to the multiplicity, \u03BC, of the pair {x, y} in the original SQS(v).")
    print("So, the multiplicity in the new system is simply \u03BC.")
    print("-" * 20)
    
    # Part (c): Must there exist pairs with multiplicity exactly v?
    v_c = 4
    # In a common doubling construction, vertical pairs {(x,0),(x,1)} have multiplicity v-1.
    mu_prime_vertical = v_c - 1
    
    # The multiplicity of a horizontal pair {(x,y)} in an SQS(v) is at most (v-2)/2.
    max_mu_horizontal = (v_c - 2) / 2
    
    print("(c) No")
    print("There is no requirement for a pair of multiplicity exactly v to exist.")
    print("Based on analysis of standard doubling constructions:")
    print(f" - For v = {v_c}, horizontal pairs {(x,0),(y,0)} have multiplicity \u03BC({x,y}), which is at most ({v_c}-2)/2 = {max_mu_horizontal:.0f}. This is less than {v_c}.")
    print(f" - Vertical pairs {(x,0),(x,1)} typically have multiplicity v-1. For v={v_c}, this is {mu_prime_vertical}, not {v_c}.")
    print("Since no class of pairs consistently results in multiplicity v, such a pair is not guaranteed to exist.")

solve_sqs_doubling_question()
print("\n<<< (a) False; (b) \u03BC; (c) No >>>")