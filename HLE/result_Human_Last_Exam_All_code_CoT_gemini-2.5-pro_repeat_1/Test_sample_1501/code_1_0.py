def solve_sqs_doubling_problem():
    """
    Solves the questions about the SQS doubling construction based on t-design properties.
    """

    # We use a valid example value for v. SQS(v) exists for v >= 4 and v congruent to 2 or 4 (mod 6).
    # Let's choose v = 10.
    v = 10

    # (a) True or False: Each element of Q x {0, 1} is contained in exactly v - 1 ND-pairs.
    # An ND-pair is of the form {(x, i), (y, i)}.
    # For a point (x, i), the other point y can be any of the v-1 points in Q other than x.
    # So, there are v-1 such pairs.
    answer_a = "True"
    print(f"(a) {answer_a}")
    print(f"For any point (x,i), there are (v-1) = {v-1} other points (y,i) to form an ND-pair.\n")


    # (b) What is the multiplicity of an ND-pair {(x, 0), (y, 0)}?
    # In an SQS(v), the multiplicity of any pair {x, y} is constant, mu = (v-2)/2.
    mu = (v - 2) / 2
    
    # In the resulting SQS(2v), the multiplicity of ANY pair is also constant, lambda_prime = (2v-2)/2 = v-1.
    lambda_prime = v - 1

    # The question asks for lambda_prime in terms of mu.
    # Since v = 2*mu + 2, then lambda_prime = v - 1 = (2*mu + 2) - 1 = 2*mu + 1.
    expression_b = "2*mu + 1"
    
    # Let's demonstrate with our v=10 example
    # mu for SQS(10) is (10-2)/2 = 4
    # lambda_prime for SQS(20) is 20-1 = 19? No, it's v-1 = 10-1 = 9
    lambda_prime_calc = v - 1
    
    # Check the expression: 2*mu + 1
    expr_val = 2 * mu + 1
    
    print(f"(b) The multiplicity of any pair in SQS(v) is mu = (v-2)/2.")
    print(f"For v={v}, mu = ({v}-2)/2 = {int(mu)}.")
    print(f"The multiplicity of any pair in SQS(2v) is lambda' = (2v-2)/2 = v-1.")
    print(f"For v={v}, lambda' = {v}-1 = {lambda_prime_calc}.")
    print(f"To express lambda' in terms of mu, we use v = 2*mu + 2.")
    print(f"So, lambda' = v - 1 = (2*mu + 2) - 1 = 2*mu + 1.")
    print(f"Checking with our numbers: {lambda_prime_calc} = 2*({int(mu)}) + 1 = {int(expr_val)}. This is correct.")
    print(f"The expression is: {expression_b}\n")
    

    # (c) Must there exist ND-pairs with multiplicity exactly v?
    # The multiplicity is v-1, which is not equal to v.
    answer_c = "No"
    print(f"(c) {answer_c}")
    print(f"The multiplicity of all ND-pairs is v-1 = {lambda_prime_calc}, which is not equal to v = {v}.\n")
    
    # Final answer format
    final_answer_string = f"<<<(a) {answer_a}; (b) {expression_b}; (c) {answer_c}>>>"
    print(final_answer_string)

solve_sqs_doubling_problem()