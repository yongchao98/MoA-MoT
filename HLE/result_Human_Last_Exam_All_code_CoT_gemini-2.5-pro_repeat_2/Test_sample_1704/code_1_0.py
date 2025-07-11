import math

def find_t4():
    """
    Finds T(4), the sum of all elements in S(4).
    S(n) is the set of all numbers that can be expressed as a sum of n
    distinct positive integers whose reciprocals sum to exactly 1.
    """
    solutions = set()
    n = 4

    # We are looking for x1, x2, x3, x4 such that:
    # 1/x1 + 1/x2 + 1/x3 + 1/x4 = 1
    # with 1 < x1 < x2 < x3 < x4

    # From 1 = sum < 4/x1, we get x1 < 4. So, x1 can be 2 or 3.
    for x1 in range(2, 4):
        # From 1 - 1/x1 = 1/x2 + 1/x3 + 1/x4 < 3/x2, we get x2 < 3*x1/(x1-1).
        x2_upper_bound = math.ceil(3 * x1 / (x1 - 1))
        for x2 in range(x1 + 1, x2_upper_bound):
            # From 1 - 1/x1 - 1/x2 = 1/x3 + 1/x4 < 2/x3, we find the upper bound for x3.
            # To avoid floating point errors, we use integer arithmetic.
            # 1/x3 + 1/x4 = (x1*x2 - x2 - x1) / (x1*x2)
            num_rem_2 = x1 * x2 - x2 - x1
            den_rem_2 = x1 * x2
            if num_rem_2 <= 0:
                continue
            x3_upper_bound = math.ceil(2 * den_rem_2 / num_rem_2)
            for x3 in range(x2 + 1, x3_upper_bound):
                # Calculate the required x4 from 1/x4 = 1 - 1/x1 - 1/x2 - 1/x3
                num_rem_3 = x1 * x2 * x3 - x2 * x3 - x1 * x3 - x1 * x2
                den_rem_3 = x1 * x2 * x3
                
                # Check if x4 is a valid integer solution
                if num_rem_3 > 0 and den_rem_3 % num_rem_3 == 0:
                    x4 = den_rem_3 // num_rem_3
                    if x4 > x3:
                        solutions.add(tuple(sorted((x1, x2, x3, x4))))

    print("The sets of 4 distinct positive integers whose reciprocals sum to 1 are:")
    s4_sums = set()
    for sol in sorted(list(solutions)):
        current_sum = sum(sol)
        s4_sums.add(current_sum)
        print(f"1/{sol[0]} + 1/{sol[1]} + 1/{sol[2]} + 1/{sol[3]} = 1, Sum = {sol[0]} + {sol[1]} + {sol[2]} + {sol[3]} = {current_sum}")
    
    t4 = sum(s4_sums)
    print("\nThe set S(4) of unique sums is:", sorted(list(s4_sums)))
    print(f"\nThe sum of all elements in S(4) is T(4) = {' + '.join(map(str, sorted(list(s4_sums))))} = {t4}")

find_t4()
<<<208>>>