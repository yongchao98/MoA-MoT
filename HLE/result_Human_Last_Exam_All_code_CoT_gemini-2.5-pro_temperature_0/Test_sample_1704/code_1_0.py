import math

def find_t4():
    """
    This function finds T(4), the sum of all elements in S(4).
    S(n) is the set of all numbers that can be expressed as a sum of n distinct
    positive integers whose reciprocals sum to exactly 1.
    """
    solutions = []
    sums = []

    # From mathematical analysis, we know x1 must be 2.
    x1 = 2

    # We also know that for x1=2, x2 must be in {3, 4, 5}.
    # We will iterate through these possibilities to find all solutions.
    # The equation is: 1/x2 + 1/x3 + 1/x4 = 1/2, with 2 < x2 < x3 < x4.

    # Loop for x2
    for x2 in range(x1 + 1, 6):
        # For a given x2, the remaining equation is 1/x3 + 1/x4 = 1/2 - 1/x2.
        # Let's represent the right-hand side as a fraction num/den.
        # 1/2 - 1/x2 = (x2 - 2) / (2 * x2)
        num = x2 - 2
        den = 2 * x2

        # We need to find x3 and x4 such that 1/x3 + 1/x4 = num/den.
        # We know x2 < x3 < x4.
        # From 1/x3 < num/den, we get x3 > den/num.
        # From 1/x3 + 1/x4 < 2/x3, we get num/den < 2/x3, so x3 < 2*den/num.
        
        x3_lower_bound = math.ceil(den / num)
        x3_upper_bound = math.floor(2 * den / num)

        for x3 in range(max(x2 + 1, x3_lower_bound), x3_upper_bound + 1):
            # Now we solve for x4: 1/x4 = num/den - 1/x3
            # 1/x4 = (num*x3 - den) / (den*x3)
            
            # The denominator for x4 must be positive
            if (num * x3 - den) <= 0:
                continue

            # Check if x4 is an integer
            if (den * x3) % (num * x3 - den) == 0:
                x4 = (den * x3) // (num * x3 - den)
                # The final check is that x4 must be distinct and larger than x3
                if x4 > x3:
                    solution_set = (x1, x2, x3, x4)
                    solutions.append(solution_set)
                    sums.append(sum(solution_set))

    # --- Output the results ---
    print("This script finds T(4), the sum of all elements in S(4).")
    print("An element of S(4) is the sum of four distinct positive integers (x1, x2, x3, x4)")
    print("whose reciprocals sum to 1 (i.e., 1/x1 + 1/x2 + 1/x3 + 1/x4 = 1).\n")
    
    print("The unique sets of integers {x1, x2, x3, x4} found are:")
    for sol in sorted(solutions):
        print(f"- {sol}, sum = {sum(sol)}")

    print("\nThe set S(4) contains the following sums:")
    print(f"S(4) = {set(sums)}\n")

    total_sum = sum(sums)
    # Sort sums in descending order for the final print statement
    sum_str = " + ".join(map(str, sorted(sums, reverse=True)))
    
    print("T(4) is the sum of all elements in S(4):")
    print(f"T(4) = {sum_str} = {total_sum}")

if __name__ == '__main__':
    find_t4()
<<<208>>>