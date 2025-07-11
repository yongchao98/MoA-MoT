import sys

def solve_archery_puzzle():
    """
    This function calculates the number of possible values for the bull's eye score
    based on the given constraints in the archery challenge.
    """
    
    # A set is used to store the unique possible values for the bull's eye score (sb).
    possible_sb_values = set()

    # The problem defines four scoring zones with strictly increasing scores: s1 < s2 < s3 < sb.
    # Scores are positive multiples of 5. Let s_i = 5 * k_i.
    # So we have integers k1 < k2 < k3 < k4, which we will call a, b, c, d.
    # This implies that d must be at least b+2, and b must be at least a+1.

    # Anna's score: 3*s1 + s2 + sb = 125.
    # In terms of our integer variables: 3*(5a) + (5b) + (5d) = 125
    # Dividing by 5 gives the core equation: 3a + b + d = 25

    # Other scores give us bounds:
    # 5*s1 <= 185 => s1 <= 37 => 5a <= 37 => a <= 7
    # 5*sb >= 230 => sb >= 46 => 5d >= 46 => d >= 9.2 => d >= 10 (since d is an integer)

    # We iterate through possible values for 'a' and 'b' to find solutions.
    
    # Loop for 'a' (scaled outer ring score). From the bound, 1 <= a <= 7.
    for a in range(1, 8):
        # Loop for 'b' (scaled second ring score).
        # We know b must be greater than a. A simple upper bound for b can be found:
        # Since d is at least 10, 3a + b + 10 <= 25 => b <= 15 - 3a.
        # So we can set the loop limit dynamically, but a fixed upper bound also works if we check all conditions.
        for b in range(a + 1, 25): 
            # Calculate 'd' (scaled bull's eye score) from Anna's equation
            d = 25 - 3 * a - b
            
            # Now, check if this set of (a, b, d) satisfies all the constraints.
            
            # Constraint 1: d must be >= 10.
            # Constraint 2: An intermediate score 's3' must exist, so d >= b + 2.
            if d >= 10 and d >= b + 2:
                # If all conditions are met, this is a valid scenario.
                # Calculate the actual bull's eye score and add it to our set.
                sb = 5 * d
                possible_sb_values.add(sb)

    # Sort the results for a clear presentation.
    sorted_values = sorted(list(possible_sb_values), reverse=True)

    print("Based on the given information, here are the possible score values for the bull's eye:")
    # The prompt requires printing each number. We will print each possible sb value.
    for value in sorted_values:
        print(value)

    print(f"\nThe exact number of possible values for the score of the bull's eye is {len(possible_sb_values)}.")

# Execute the function to solve the puzzle.
solve_archery_puzzle()

sys.stdout.flush()