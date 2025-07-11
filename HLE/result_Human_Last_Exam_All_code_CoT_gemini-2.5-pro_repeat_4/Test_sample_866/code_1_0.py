import math

def solve_crease_pattern():
    """
    Calculates the number of valid flat-foldable assignments for a given
    single-vertex crease pattern.
    """
    # The input pattern as a Python list.
    pattern = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']

    # 1. Parse the input into angles and creases.
    angles = [item for i, item in enumerate(pattern) if i % 2 == 0]
    creases = [item for i, item in enumerate(pattern) if i % 2 != 0]

    # 2. Check the most fundamental condition for flat-folding: the number of creases.
    # For a single vertex to be flat-foldable, it must have an even number of creases.
    num_creases = len(creases)
    
    print(f"The given pattern has {num_creases} creases.")

    # 3. If the number of creases is odd, no solution exists.
    if num_creases % 2 != 0:
        print(f"Checking if the number of creases is even: {num_creases} % 2 = {num_creases % 2}")
        print("The result is not 0, so the number of creases is odd.")
        print("A single vertex with an odd number of creases cannot be folded flat.")
        count = 0
    else:
        # This part of the code would run if the number of creases were even.
        # It is included for completeness but will not be executed for the given input.
        
        # Check Kawasaki's Theorem: sum of alternate angles must be 180 degrees.
        # Note: round() is used to handle potential floating-point inaccuracies.
        sum_alt1 = sum(angles[i] for i in range(0, num_creases, 2))
        sum_alt2 = sum(angles[i] for i in range(1, num_creases, 2))

        if round(sum_alt1) != 180 or round(sum_alt2) != 180:
            print("Kawasaki's Theorem is not satisfied.")
            count = 0
        else:
            # Check Maekawa's Theorem: |#M - #V| = 2
            num_M_known = creases.count('M')
            num_V_known = creases.count('V')
            num_unknown = creases.count('?')
            
            count = 0
            # Iterate through all possibilities of assigning 'M' to the '?' creases.
            for k in range(num_unknown + 1):
                # k is the number of '?' assigned to 'M'
                total_M = num_M_known + k
                total_V = num_V_known + (num_unknown - k)
                
                if abs(total_M - total_V) == 2:
                    # If Maekawa's is satisfied, count the number of ways to choose k '?'s.
                    # This is C(num_unknown, k).
                    combinations = math.comb(num_unknown, k)
                    count += combinations

    print(f"\nTotal number of different flat-foldable assignments: {count}")

solve_crease_pattern()