import math

def find_t4():
    """
    Finds T(4), the sum of all elements in S(4), where S(n) is the set of
    sums of n distinct positive integers whose reciprocals sum to 1.
    This is done by finding all sets {x1, x2, x3, x4} satisfying
    1/x1 + 1/x2 + 1/x3 + 1/x4 = 1.
    """
    solutions = set()
    epsilon = 1e-9

    # Assume 1 < x1 < x2 < x3 < x4 to find unique sets of integers.
    # From the inequality 1 < 4/x1, we deduce x1 < 4. So, x1 is 2 or 3.
    for x1 in range(2, 4):
        # From 1 - 1/x1 < 3/x2, we find an upper limit for x2.
        x2_limit = 3 * x1 / (x1 - 1)
        for x2 in range(x1 + 1, math.ceil(x2_limit)):
            
            rem_after_x2 = 1.0 - 1.0/x1 - 1.0/x2
            if rem_after_x2 <= epsilon:
                continue

            # We derive lower and upper bounds for x3 to narrow the search space.
            x3_lower_bound = 1.0 / rem_after_x2
            x3_limit = 2.0 / rem_after_x2
            
            start_x3 = max(x2 + 1, math.floor(x3_lower_bound) + 1)
            end_x3 = math.ceil(x3_limit)
            
            for x3 in range(start_x3, end_x3):
                rem_after_x3 = 1.0 - 1.0/x1 - 1.0/x2 - 1.0/x3
                
                if rem_after_x3 <= epsilon:
                    continue

                # Calculate the potential value of x4.
                x4_float = 1.0 / rem_after_x3
                x4_int = int(round(x4_float))
                
                # Check if x4 is a distinct integer.
                if abs(x4_int - x4_float) < epsilon and x4_int > x3:
                    solution_tuple = tuple(sorted((x1, x2, x3, x4_int)))
                    solutions.add(solution_tuple)

    # S(4) is the set of sums of these integer solutions.
    s4 = sorted(list({sum(s) for s in solutions}))

    # T(4) is the sum of all elements in S(4).
    t4 = sum(s4)

    # Print the final equation as requested.
    sum_str = " + ".join(map(str, s4))
    print(f"{sum_str} = {t4}")

find_t4()