def solve_sqs_problem():
    """
    Solves the theoretical problem about the SQS doubling construction based on provided definitions.
    The solution relies on logical deduction from the properties of the designs rather than a specific
    implementation of the construction.
    """

    # Step 1: Establish properties from definitions.
    # Initial SQS(v): The multiplicity mu of any pair is given by the parameter
    # of the 2-(v, 2, v-1) design, so mu = v - 1.

    # Final SQS(2v): The ND-pairs form a 2-(2v, 2, 2v-1) design.
    # Let n = 2v, k = 2, lambda_prime = 2v - 1.

    # Step 2: Solve part (a).
    # The number of ND-pairs containing any single element is the replication number 'r_prime'.
    # r_prime = lambda_prime * (n - 1) / (k - 1)
    # r_prime = (2v - 1) * (2v - 1) / (2 - 1) = (2v - 1).
    # The question asks if this value is v - 1.
    # For v >= 4, (2v - 1) is not equal to (v - 1).
    # Therefore, the statement is False.
    answer_a = "False"

    # Step 3: Solve part (b).
    # The multiplicity of any ND-pair in the new system is lambda_prime = 2v - 1.
    # We need to express this in terms of mu.
    # From the initial system, we have v = mu + 1.
    # Substitute this into the expression for the new multiplicity:
    # Multiplicity = 2 * (mu + 1) - 1
    #             = 2*mu + 2 - 1
    #             = 2*mu + 1
    # The numbers in this final equation are 2 and 1.
    answer_b_expression = "2*mu + 1"

    # Step 4: Solve part (c).
    # The multiplicity of all ND-pairs is 2v - 1.
    # We check if it's possible for this to equal v.
    # 2v - 1 = v  =>  v = 1.
    # This contradicts the given condition v >= 4.
    # Therefore, no ND-pair can have multiplicity v.
    answer_c = "No"

    # Step 5: Print the final answer in the specified format.
    print(f"(a) [{answer_a}]; (b) [{answer_b_expression}]; (c) [{answer_c}].")

solve_sqs_problem()