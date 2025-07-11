import math
from fractions import Fraction

def find_t4():
    """
    Finds all sets of 4 distinct positive integers whose reciprocals sum to 1,
    calculates the sum of each set's elements (S(4)), and then finds the sum
    of all those unique sums (T(4)).
    """
    solutions = []
    s_values = set()

    # We are looking for x1, x2, x3, x4 such that:
    # 1/x1 + 1/x2 + 1/x3 + 1/x4 = 1
    # with 0 < x1 < x2 < x3 < x4

    # Bound for x1: 1 = 1/x1 + 1/x2 + 1/x3 + 1/x4 < 4/x1 => x1 < 4.
    # x1 cannot be 1, so x1 is in {2, 3}.
    for x1 in range(2, 4):
        # Remaining sum needed after choosing x1
        rem1 = Fraction(1) - Fraction(1, x1)

        # Bound for x2: rem1 = 1/x2 + 1/x3 + 1/x4 < 3/x2 => x2 < 3 / rem1
        x2_upper_bound = 3 / rem1
        for x2 in range(x1 + 1, math.ceil(x2_upper_bound)):
            # Remaining sum needed after choosing x2
            rem2 = rem1 - Fraction(1, x2)
            if rem2 <= 0:
                continue

            # Bound for x3: rem2 = 1/x3 + 1/x4 < 2/x3 => x3 < 2 / rem2
            # Also, rem2 > 1/x3 => x3 > 1 / rem2
            x3_upper_bound = 2 / rem2
            x3_lower_bound = 1 / rem2
            
            start_x3 = max(x2 + 1, math.floor(x3_lower_bound) + 1)

            for x3 in range(start_x3, math.ceil(x3_upper_bound)):
                # Remaining sum for 1/x4
                rem3 = rem2 - Fraction(1, x3)

                # Check if rem3 can be a reciprocal of an integer
                if rem3 > 0 and rem3.numerator == 1:
                    x4 = rem3.denominator
                    # Ensure x4 is distinct and ordered
                    if x4 > x3:
                        solution_set = (x1, x2, x3, x4)
                        solutions.append(solution_set)
                        s_values.add(sum(solution_set))

    print("Found the following solutions {x1, x2, x3, x4} and their sums:")
    # Sort solutions for consistent output
    solutions.sort()
    for sol in solutions:
        s = sum(sol)
        print(f"1/{sol[0]} + 1/{sol[1]} + 1/{sol[2]} + 1/{sol[3]} = 1  =>  Sum = {sol[0]} + {sol[1]} + {sol[2]} + {sol[3]} = {s}")

    print("\nThe set of sums, S(4), is:", sorted(list(s_values)))

    total_sum_t4 = sum(s_values)
    
    # Format the final equation as requested
    sorted_s = sorted(list(s_values))
    equation_str = " + ".join(map(str, sorted_s))
    
    print("\nThe sum of elements in S(4), T(4), is:")
    print(f"T(4) = {equation_str} = {total_sum_t4}")
    
    print(f"\n<<<T(4) = {total_sum_t4}>>>")

if __name__ == "__main__":
    find_t4()