import math

def find_t4():
    """
    Finds all sets of 4 distinct positive integers whose reciprocals sum to 1,
    calculates the sum of each set (elements of S(4)), and then computes
    the sum of these sums (T(4)).
    """
    solutions = []
    s4_elements = set()

    # We assume x1 < x2 < x3 < x4.
    # From 1 = 1/x1 + 1/x2 + 1/x3 + 1/x4 < 4/x1, we get x1 < 4.
    # If x1 = 1, 1/x1 = 1, which is impossible for distinct positive integers.
    # So, 2 <= x1 <= 3.
    for x1 in range(2, 4):
        # From 1 - 1/x1 = 1/x2 + 1/x3 + 1/x4 < 3/x2, we get x2 < 3 / (1 - 1/x1).
        upper_bound_x2 = 3 / (1 - 1/x1)
        for x2 in range(x1 + 1, math.ceil(upper_bound_x2)):
            # From 1 - 1/x1 - 1/x2 = 1/x3 + 1/x4 < 2/x3, we get x3 < 2 / (1 - 1/x1 - 1/x2).
            # Also, 1/x3 < 1 - 1/x1 - 1/x2, so x3 > 1 / (1-1/x1-1/x2)
            rem_sum_2 = 1 - 1/x1 - 1/x2
            if rem_sum_2 <= 0:
                continue
            
            lower_bound_x3 = 1 / rem_sum_2
            upper_bound_x3 = 2 / rem_sum_2
            
            for x3 in range(max(x2 + 1, math.ceil(lower_bound_x3)), math.ceil(upper_bound_x3)):
                # Calculate x4 from 1/x4 = 1 - 1/x1 - 1/x2 - 1/x3
                # To avoid floating point issues, we work with numerators and denominators.
                # 1/x4 = (x1*x2*x3 - x2*x3 - x1*x3 - x1*x2) / (x1*x2*x3)
                num_x4 = x1 * x2 * x3
                den_x4 = num_x4 - (x2 * x3) - (x1 * x3) - (x1 * x2)

                if den_x4 > 0 and num_x4 % den_x4 == 0:
                    x4 = num_x4 // den_x4
                    if x4 > x3:
                        solution_set = tuple(sorted((x1, x2, x3, x4)))
                        current_sum = sum(solution_set)
                        if current_sum not in s4_elements:
                            solutions.append(solution_set)
                            s4_elements.add(current_sum)
    
    print("Found the following sets {x1, x2, x3, x4} and their sums (s):")
    for s_set in sorted(solutions, key=lambda x: sum(x)):
        s = sum(s_set)
        print(f"Set = {s_set}, Sum = {s}")

    print("\nThe set S(4) is the set of these sums:")
    sorted_s4 = sorted(list(s4_elements))
    print(f"S(4) = {set(sorted_s4)}")
    
    t4 = sum(sorted_s4)
    sum_equation = " + ".join(map(str, sorted_s4))
    
    print("\nT(4) is the sum of the elements in S(4):")
    print(f"T(4) = {sum_equation} = {t4}")

if __name__ == "__main__":
    find_t4()