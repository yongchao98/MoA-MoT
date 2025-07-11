def solve_sqs_doubling_problem():
    """
    Solves and explains the answers to the nested SQS doubling construction problem.
    """
    print("This solution is based on the following interpretations of the terms:")
    print("- A nested SQS(v) for v even has neighbor pairs forming a 1-factorization of K_v.")
    print("- An 'ND-pair' is a neighbor pair in the resulting SQS(2v) of the form {(x, 0), (y, 1)}.")
    print("- 'Multiplicity' of a pair is the number of blocks it is contained in.\n")

    # --- Part (a) ---
    print("--- Part (a) ---")
    print("Question: True or False: In the doubling construction, each element of Q x {0, 1} is contained in exactly v - 1 ND-pairs.")
    print("Answer: False")
    print("\nReasoning:")
    print("The neighbor pairs in the resulting nested SQS(2v) form a 1-factorization of the complete graph K_2v.")
    print("This means every point is connected to every other point by exactly one neighbor pair.")
    print("Consider an element p = (x, 0). It must form neighbor pairs with all 2v-1 other points.")
    print("These other points consist of:")
    print("  - v-1 points of the form (y, 0) where y != x.")
    print("  - v points of the form (y, 1).")
    print("An ND-pair containing p must be of the form {(x, 0), (y, 1)}. There are exactly 'v' such pairs.")
    print("Since v != v-1, the statement is false.\n")

    # --- Part (b) ---
    print("--- Part (b) ---")
    print("Question: What is the multiplicity of an ND-pair {(x, 0), (y, 0)} in the resulting nested SQS(2v) if the pair {x, y} had multiplicity mu in the original SQS(v)?")
    print("Answer: The expression is mu + v/2")
    print("\nReasoning:")
    print("The question contains a typo, as {(x, 0), (y, 0)} is not an ND-pair. We assume it means 'pair'.")
    print("The multiplicity of any pair in an SQS(2v) is (2v-2)/2 = v-1.")
    print("The multiplicity of the pair {x, y} in the original SQS(v) is mu = (v-2)/2.")
    print("A standard doubling construction embeds the original SQS(v) into the new structure. This accounts for 'mu' of the total multiplicity for the pair {(x, 0), (y, 0)}.")
    print("The remaining multiplicity must come from new blocks created by the construction. This amount is (v-1) - mu.")
    print("(v-1) - mu = (v-1) - (v-2)/2 = (2v-2 - v+2)/2 = v/2.")
    print("Thus, the total multiplicity can be expressed as the sum of these parts: mu + v/2.\n")
    
    print("Example for v=4 (the smallest possible value):")
    v = 4
    mu = (v - 2) / 2
    v_div_2 = v / 2
    total_multiplicity = mu + v_div_2
    print(f"For v = {v}, mu = ({v}-2)/2 = {int(mu)}.")
    print(f"The second term is v/2 = {v}/2 = {int(v_div_2)}.")
    print("The final equation is:")
    print(f"{int(mu)} + {int(v_div_2)} = {int(total_multiplicity)}")
    print(f"This matches the expected multiplicity in an SQS(8), which is (8-2)/2 = 3.\n")

    # --- Part (c) ---
    print("--- Part (c) ---")
    print("Question: Must there exist ND-pairs with multiplicity exactly v in the constructed nested SQS(2v)?")
    print("Answer: No")
    print("\nReasoning:")
    print("An ND-pair is a pair of points, e.g., {(x, 0), (y, 1)}.")
    print("The multiplicity of any pair of points in a valid SQS(2v) is strictly defined as (2v-2)/2 = v-1.")
    print("The multiplicity cannot be v, as this would violate the definition of a Steiner Quadruple System.")
    print("Therefore, no pair, ND or otherwise, can have multiplicity v.\n")

    # --- Final Answer Summary ---
    print("--- Summary ---")
    print("(a) [False]; (b) [mu + v/2]; (c) [No]")


if __name__ == '__main__':
    solve_sqs_doubling_problem()