def find_t4():
    """
    Finds T(4), the sum of all elements in S(4). S(4) is the set of all numbers 
    that can be expressed as a sum of 4 distinct positive integers whose 
    reciprocals sum to exactly 1.
    """
    solutions = []
    # We are looking for distinct positive integers a, b, c, d such that:
    # 1/a + 1/b + 1/c + 1/d = 1
    # We assume 1 < a < b < c < d to find unique sets.
    
    # Based on analysis, a can only be 2 or 3.
    # Loop for a
    for a in range(2, 4):
        # Bound for b: 3/b > 1 - 1/a => b < 3a/(a-1)
        b_upper_bound = (3 * a) // (a - 1) + 1
        for b in range(a + 1, b_upper_bound):
            # Bound for c: 2/c > 1 - 1/a - 1/b => c < 2ab/(ab-a-b)
            # Check if denominator is positive
            if a * b - a - b <= 0:
                continue
            c_upper_bound = (2 * a * b) // (a * b - a - b) + 1
            for c in range(b + 1, c_upper_bound):
                # Calculate d from 1/d = 1 - 1/a - 1/b - 1/c
                # d = abc / (abc - bc - ac - ab)
                numerator_d_inv = a * b * c - b * c - a * c - a * b
                if numerator_d_inv > 0:
                    denominator_d_inv = a * b * c
                    # Check if d is an integer
                    if denominator_d_inv % numerator_d_inv == 0:
                        d = denominator_d_inv // numerator_d_inv
                        # Check if d is distinct and follows the order
                        if d > c:
                            solutions.append((a, b, c, d))

    print("Finding all sets of 4 distinct positive integers {a, b, c, d} such that 1/a + 1/b + 1/c + 1/d = 1.")
    print("Assuming a < b < c < d, we systematically search for solutions.")

    # Sort solutions based on their sum for clearer presentation
    solutions.sort(key=sum)

    # Calculate the set of sums S(4)
    s4 = set()
    print("\nThe sets of integers found and their sums are:")
    for s in solutions:
        current_sum = sum(s)
        s4.add(current_sum)
        print(f"Set: {s}, Sum: {s[0]} + {s[1]} + {s[2]} + {s[3]} = {current_sum}")

    # Sort the unique sums for the final equation
    s4_list = sorted(list(s4))

    print("\nThe set S(4) of unique sums is:")
    print(f"S(4) = {s4_list}")
    
    # Calculate T(4)
    total_sum = sum(s4_list)
    
    print("\nT(4) is the sum of the elements in S(4):")
    # Build the equation string showing each number
    equation_str = " + ".join(map(str, s4_list))
    print(f"T(4) = {equation_str} = {total_sum}")

    # Return the final numerical answer as requested by the format.
    # The print statements above provide the detailed explanation.
    return total_sum

if __name__ == '__main__':
    # This block will be executed when the script is run.
    T4 = find_t4()
    # The final answer is wrapped according to the user's request.
    # print(f"<<<{T4}>>>")
    # For this interactive session, I will just let the function print everything.
    # The final answer will be manually placed at the very end.

find_t4()