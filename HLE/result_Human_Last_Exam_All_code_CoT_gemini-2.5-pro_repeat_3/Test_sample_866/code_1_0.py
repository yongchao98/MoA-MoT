import math

def solve_flat_folding():
    """
    Calculates the number of valid crease assignments for a given single-vertex
    pattern to be flat-foldable.
    """
    # The input pattern: [angle, crease, angle, crease, ...]
    pattern = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']

    # Step 1: Extract the crease assignments from the pattern.
    # Creases are at odd-numbered indices (1, 3, 5, ...).
    creases = [item for i, item in enumerate(pattern) if i % 2 != 0]

    # Step 2: Count the total number of creases.
    n = len(creases)
    print(f"The given pattern has a total of {n} creases.")

    # Step 3: Apply Maekawa's Theorem and its consequences.
    # Maekawa's Theorem states that for a flat-foldable vertex, the number of
    # mountain folds (M) and valley folds (V) must satisfy the equation:
    # abs(M - V) = 2
    # The total number of creases is n = M + V.
    #
    # A key consequence is that the total number of creases (n) must be even.
    # Proof: The sum (M+V) and difference (M-V) of two integers must have
    # the same parity. Since |M-V|=2 is even, n = M+V must also be even.

    if n % 2 != 0:
        print(f"The total number of creases is {n}, which is an odd number.")
        print("According to Maekawa's theorem, a flat-foldable vertex must have an even number of creases.")
        print("Therefore, no combination of assignments for '?' can satisfy the condition |M - V| = 2.")
        final_answer = 0
    else:
        # This part of the code would run if n were even.
        known_M = creases.count('M')
        known_V = creases.count('V')
        unknown_count = creases.count('?')
        
        total_assignments = 0
        
        # Iterate through all ways to assign 'M' to the '?' creases.
        for m_q in range(unknown_count + 1):
            v_q = unknown_count - m_q
            
            total_M = known_M + m_q
            total_V = known_V + v_q
            
            # Check if Maekawa's theorem holds for this assignment.
            if abs(total_M - total_V) == 2:
                # If it holds, count the number of ways to choose m_q '?'s to be 'M'.
                num_ways = math.comb(unknown_count, m_q)
                print(f"Found a valid assignment: {total_M} Mountains and {total_V} Valleys.")
                print(f"Equation: abs({total_M} - {total_V}) = {abs(total_M - total_V)}")
                print(f"Number of ways for this assignment: {num_ways}")
                total_assignments += num_ways
        final_answer = total_assignments

    print(f"\nTotal number of different flat-foldable assignments: {final_answer}")

solve_flat_folding()