def find_t4():
    """
    This function finds T(4), the sum of all numbers that can be expressed as a sum of 4
    distinct positive integers whose reciprocals sum to exactly 1.
    """
    n = 4
    solutions = []
    sums = set()

    # We are looking for distinct positive integers x1, x2, x3, x4 such that:
    # 1/x1 + 1/x2 + 1/x3 + 1/x4 = 1
    # Assume x1 < x2 < x3 < x4.

    # Bound for x1: 1 = 1/x1 + ... < n/x1 => x1 < n. Since x1 > 1, x1 can be 2 or 3.
    for x1 in range(2, n):
        # Bound for x2: 1 - 1/x1 < (n-1)/x2 => x2 < (n-1)*x1/(x1-1)
        x2_ub_num = (n - 1) * x1
        x2_ub_den = x1 - 1
        x2_upper_limit = x2_ub_num // x2_ub_den

        for x2 in range(x1 + 1, x2_upper_limit + 1):
            # Ensure strict inequality for the bound
            if x2 * x2_ub_den >= x2_ub_num:
                continue

            # Remainder for x3, x4: 1/x3 + 1/x4 = 1 - 1/x1 - 1/x2
            # We use integer arithmetic to avoid float precision issues.
            # rem2 = (x2*(x1-1) - x1) / (x1*x2)
            num2 = x2 * (x1 - 1) - x1
            den2 = x1 * x2
            if num2 <= 0:
                continue

            # Bound for x3: num2/den2 < (n-2)/x3 => x3 < (n-2)*den2/num2
            x3_ub_num = (n - 2) * den2
            x3_ub_den = num2
            x3_upper_limit = x3_ub_num // x3_ub_den
            
            # Lower bound for x3: 1/x3 < num2/den2 => x3 > den2/num2
            x3_lower_limit = den2 // num2 + 1
            
            x3_start = max(x2 + 1, x3_lower_limit)

            for x3 in range(x3_start, x3_upper_limit + 1):
                # Ensure strict inequality for the bound
                if x3 * x3_ub_den >= x3_ub_num:
                    continue

                # Remainder for x4: 1/x4 = num2/den2 - 1/x3
                # rem3 = (x3*num2 - den2) / (den2*x3)
                num3 = x3 * num2 - den2
                den3 = den2 * x3
                if num3 <= 0:
                    continue

                # Check if x4 is an integer
                if den3 % num3 == 0:
                    x4 = den3 // num3
                    # Check if x4 is distinct and maintains the order
                    if x4 > x3:
                        solution = (x1, x2, x3, x4)
                        solutions.append(solution)
                        sums.add(sum(solution))

    print("Found the following sets of 4 distinct positive integers whose reciprocals sum to 1:")
    # Sort solutions for consistent output
    solutions.sort()
    for sol in solutions:
        s = sum(sol)
        print(f"- Set: {sol}, Sum: {s}")
        print(f"  Check: 1/{sol[0]} + 1/{sol[1]} + 1/{sol[2]} + 1/{sol[3]} = 1")
    
    sorted_sums = sorted(list(sums))
    total_sum = sum(sorted_sums)

    print("\nThe set S(4) is the set of these unique sums:")
    print(f"S(4) = {set(sorted_sums)}")
    
    print("\nT(4) is the sum of the elements in S(4). The final equation is:")
    equation = " + ".join(map(str, sorted_sums))
    print(f"T(4) = {equation} = {total_sum}")

if __name__ == '__main__':
    find_t4()