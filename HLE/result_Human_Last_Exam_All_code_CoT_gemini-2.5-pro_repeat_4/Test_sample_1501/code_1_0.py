def solve_combinatorial_problem():
    """
    This script analyzes the properties of a nested Steiner Quadruple System (SQS)
    created via a doubling construction and answers the specific questions posed.
    """
    
    # Introduction to the reasoning process
    print("This program solves the problem by following a logical deduction based on the properties of SQS and the specified doubling construction.")
    print("-----------------------------------------------------------------------------------------------------------------------------------")
    
    # Step 1 & 2: Define the construction and key terms
    print("\nStep 1: Understanding the Doubling Construction and 'ND-pair'")
    print("Let S = (Q, B) be an SQS(v). The point set of the new SQS(2v) is Q' = Q x {0, 1}.")
    print("The blocks B' of the SQS(2v) are formed by two types:")
    print(" - Type 1: For each block {a,b,c,d} in B, we form 8 blocks {(a,i),(b,j),(c,k),(d,l)} where i,j,k,l in {0,1} and i+j+k+l is even.")
    print(" - Type 2: For each pair {x,y} from Q, we form the block {(x,0),(y,0),(x,1),(y,1)}.")
    print("\nThe problem states each point in Q' is in v-1 'ND-pairs'. We deduce that an ND-pair is a 'horizontal' pair: {(x, a), (y, a)} for x != y.")
    print("A point (x,a) is part of v-1 such pairs, one for each y in Q other than x.")
    
    # Step 3: Analyze and answer each question
    print("\nStep 2: Answering the questions based on the definitions.")

    # (a) True or False: Each element is in exactly v-1 ND-pairs.
    answer_a = "True"
    print(f"\n(a) Question: Is it '{answer_a}' that each element is in exactly v-1 ND-pairs?")
    print(f"   Reasoning: As established, an element (x, a) forms an ND-pair with every other element (y, a) on the same 'level'. There are v-1 such elements y.")
    print(f"   Conclusion: The statement is {answer_a}.")

    # (b) What is the multiplicity of an ND-pair?
    # The pair {x, y} in the SQS(v) has multiplicity μ (mu).
    # We calculate the multiplicity of the ND-pair {(x, 0), (y, 0)} in the SQS(2v).
    # Contribution from Type 1 blocks: 2 * μ
    # Contribution from Type 2 blocks: 1
    # Total multiplicity = 2 * μ + 1
    answer_b = "2\u03BC + 1"
    print(f"\n(b) Question: What is the multiplicity of an ND-pair {(x, 0), (y, 0)} if the original pair {{x, y}} had multiplicity \u03BC?")
    print(f"   Reasoning: We count the blocks in SQS(2v) containing {(x, 0), (y, 0)}:")
    print(f"    - For each of the \u03BC blocks in SQS(v) containing {{x,y}}, two Type 1 blocks are generated. Contribution: 2 * \u03BC.")
    print(f"    - The pair {{x,y}} defines exactly one Type 2 block. Contribution: 1.")
    print(f"   Conclusion: The total multiplicity is the sum of the contributions.")
    print(f"   Final expression for multiplicity: {answer_b}")

    # (c) Must there exist ND-pairs with multiplicity exactly v?
    # In a simple SQS(v), the multiplicity μ is constant for all pairs: μ = (v-2)/2.
    # So, the multiplicity of any ND-pair in the SQS(2v) is 2 * ((v-2)/2) + 1 = v - 1.
    # Can v-1 = v? No, for v >= 4.
    answer_c = "No"
    print(f"\n(c) Question: Must there exist ND-pairs with multiplicity exactly v?")
    print(f"   Reasoning: The multiplicity of an ND-pair is 2\u03BC + 1. For any SQS(v), it is a known theorem that \u03BC must equal (v-2)/2.")
    print(f"   Therefore, the multiplicity of every ND-pair is 2*((v-2)/2) + 1 = (v-2) + 1 = v - 1.")
    print(f"   The condition 'multiplicity = v' becomes v-1 = v, which is impossible.")
    print(f"   Conclusion: The answer is {answer_c}.")

    print("-----------------------------------------------------------------------------------------------------------------------------------")
    # Consolidate and print the final answer in the required format
    final_answer_string = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
    print("\nFinal Answer:")
    print(final_answer_string)


# Execute the function to solve the problem
solve_combinatorial_problem()