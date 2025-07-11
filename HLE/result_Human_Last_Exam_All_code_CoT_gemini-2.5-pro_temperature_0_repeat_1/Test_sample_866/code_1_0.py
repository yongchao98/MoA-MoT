import math

def solve_flat_fold():
    """
    Calculates the number of valid flat-foldable assignments for a given
    partially assigned single-vertex crease pattern.
    The function prints the final calculation or result.
    """
    # The input crease pattern provided by the user.
    pattern_data = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']

    # Step 1: Parse the input into separate lists for angles and creases.
    angles = [item for item in pattern_data if isinstance(item, (int, float))]
    creases = [item for item in pattern_data if isinstance(item, str)]

    num_creases = len(creases)

    # Step 2: Check the fundamental conditions for flat-folding.

    # Condition A: The total number of creases must be even.
    # This is a requirement for both Maekawa's and Kawasaki's theorems.
    if num_creases % 2 != 0:
        # If the number of creases is odd, no flat-foldable assignment is possible.
        # The number of valid assignments is 0. There is no equation to compute.
        print(0)
        return

    # Condition B: Kawasaki's Theorem requires the sum of alternating angles to be 180°.
    # This part of the code is not reached for the given input but is included for a general solution.
    sum_odd_indexed_angles = sum(angles[i] for i in range(0, num_creases, 2))
    
    # We use math.isclose for robust floating-point comparison.
    if not math.isclose(sum_odd_indexed_angles, 180):
        # If the angles do not satisfy the condition, no assignment is valid.
        print(0)
        return

    # Step 3: If the above conditions are met, count assignments that satisfy Maekawa's Theorem.
    num_known_m = creases.count('M')
    num_known_v = creases.count('V')
    num_unknown = creases.count('?')

    # This list will store the number of ways for each valid scenario.
    combination_counts = []

    # Iterate through all possible ways to assign 'M' to the '?' creases.
    # 'k' represents the number of '?'s that are assigned as 'M'.
    for k in range(num_unknown + 1):
        total_m = num_known_m + k
        total_v = num_known_v + (num_unknown - k)

        # Check Maekawa's Theorem: |#M - #V| = 2
        if abs(total_m - total_v) == 2:
            # If the assignment is valid, calculate the number of ways to achieve it.
            # This is the binomial coefficient C(num_unknown, k).
            count = math.comb(num_unknown, k)
            combination_counts.append(count)

    total_assignments = sum(combination_counts)

    # Step 4: Print the final result, including the equation if applicable.
    if total_assignments == 0:
        # This occurs if no combination of '?' assignments satisfies Maekawa's theorem.
        print(0)
    else:
        # To meet the requirement, we show the numbers that sum up to the final answer.
        equation_str = " + ".join(map(str, combination_counts))
        print(f"{equation_str} = {total_assignments}")

# Execute the solver function.
solve_flat_fold()