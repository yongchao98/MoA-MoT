import math

def find_t4():
    """
    Finds all sets of 4 distinct positive integers whose reciprocals sum to 1,
    calculates the sum of integers for each set (S(4)),
    and then calculates the sum of those sums (T(4)).
    """
    n = 4
    solutions = []
    found_sums = set()

    # We assume x1 < x2 < x3 < x4 to find unique sets.
    # From 1 = 1/x1 + 1/x2 + 1/x3 + 1/x4 < 4/x1, we get x1 < 4.
    # Since x1 must be > 1, x1 can be 2 or 3.
    for x1 in range(2, 4):
        # From 1 - 1/x1 = 1/x2 + 1/x3 + 1/x4 < 3/x2, we get x2 < 3*x1/(x1-1).
        x2_upper_bound = math.ceil(3 * x1 / (x1 - 1))
        for x2 in range(x1 + 1, x2_upper_bound):
            # Let rem_frac = 1 - 1/x1 - 1/x2.
            # From rem_frac = 1/x3 + 1/x4 < 2/x3, we get x3 < 2/rem_frac.
            # Using integer arithmetic for rem_frac = (x1*x2 - x2 - x1) / (x1*x2)
            rem_num = x1 * x2 - x2 - x1
            rem_den = x1 * x2
            if rem_num <= 0:
                continue
            
            x3_upper_bound = math.ceil(2 * rem_den / rem_num)
            # The lower bound for x3 is x2 + 1.
            # Also x3 > 1/rem_frac = rem_den/rem_num.
            x3_lower_bound = max(x2 + 1, math.floor(rem_den / rem_num) + 1)
            
            for x3 in range(x3_lower_bound, x3_upper_bound):
                # Calculate x4 from 1/x4 = rem_frac - 1/x3
                # 1/x4 = (rem_num * x3 - rem_den) / (rem_den * x3)
                num_x4 = rem_den * x3
                den_x4 = rem_num * x3 - rem_den

                if den_x4 > 0 and num_x4 % den_x4 == 0:
                    x4 = num_x4 // den_x4
                    if x4 > x3:
                        solution_set = (x1, x2, x3, x4)
                        current_sum = sum(solution_set)
                        if current_sum not in found_sums:
                            found_sums.add(current_sum)
                            solutions.append(solution_set)

    print("Found the following sets {x1, x2, x3, x4} and their sums:")
    # Sort solutions based on their sum for clear presentation
    solutions.sort(key=sum)
    for sol in solutions:
        print(f"Set: {{{', '.join(map(str, sol))}}}, Sum = {sum(sol)}")

    print("\nThe set of sums S(4) is:")
    sum_list = sorted(list(found_sums))
    print(f"S(4) = {{{', '.join(map(str, sum_list))}}}")

    print("\nT(4) is the sum of the elements in S(4).")
    total_sum = sum(sum_list)
    equation = " + ".join(map(str, sum_list))
    print(f"T(4) = {equation} = {total_sum}")

if __name__ == '__main__':
    find_t4()