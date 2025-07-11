import math

def find_t4():
    """
    Solves for T(4), which is the sum of all elements in S(4).
    S(n) is the set of sums of n distinct positive integers whose reciprocals sum to 1.
    This function finds all sets {x1, x2, x3, x4} satisfying the condition,
    calculates their sums, and then sums those unique sums.
    """
    s4_sums = set()
    solutions = []

    print("Step 1: Find all sets of 4 distinct positive integers {x1, x2, x3, x4} such that:")
    print("1/x1 + 1/x2 + 1/x3 + 1/x4 = 1\n")

    # Assuming 1 < x1 < x2 < x3 < x4.
    # From analysis, x1 can only be 2 or 3.
    for x1 in range(2, 4):
        # Upper bound for x2 from 1 - 1/x1 < 3/x2
        x2_upper_bound = math.ceil((3 * x1) / (x1 - 1))
        for x2 in range(x1 + 1, x2_upper_bound):
            # Let R = 1 - 1/x1 - 1/x2 = N/D
            N = x1 * x2 - x1 - x2
            D = x1 * x2
            if N <= 0:
                continue

            # Upper bound for x3 from R < 2/x3 -> x3 < 2*D/N
            # Lower bound for x3 from R > 1/x3 -> x3 > D/N
            x3_lower_bound = math.floor(D / N) + 1
            x3_upper_bound = math.ceil((2 * D) / N)
            
            start_x3 = max(x2 + 1, x3_lower_bound)

            for x3 in range(start_x3, x3_upper_bound):
                # Solve for x4 from 1/x4 = R - 1/x3 = (N*x3 - D) / (D*x3)
                # which means x4 = (D*x3) / (N*x3 - D)
                num_x4 = D * x3
                den_x4 = N * x3 - D

                if den_x4 > 0 and num_x4 % den_x4 == 0:
                    x4 = num_x4 // den_x4
                    if x4 > x3:
                        # Found a valid solution set
                        solution_set = sorted([x1, x2, x3, x4])
                        solutions.append(solution_set)
                        s4_sums.add(sum(solution_set))

    print("Step 2: List all found solution sets and their sums.")
    print("-" * 50)
    sorted_solutions = sorted(solutions)
    for sol in sorted_solutions:
        s = sum(sol)
        print(f"Found solution set: {sol}")
        print(f"Equation: {sol[0]} + {sol[1]} + {sol[2]} + {sol[3]} = {s}")
        print("-" * 50)

    print("Step 3: Determine the set of unique sums, S(4).")
    sorted_sums = sorted(list(s4_sums))
    print(f"S(4) = {sorted_sums}\n")

    print("Step 4: Calculate T(4) by summing the elements of S(4).")
    total_sum = sum(s4_sums)
    sum_equation = " + ".join(map(str, sorted_sums))
    print(f"T(4) = {sum_equation} = {total_sum}")

if __name__ == '__main__':
    find_t4()