import itertools

def solve_flat_folding():
    """
    Determines the number of valid flat-foldable assignments for a given
    partially assigned crease pattern at a single vertex.
    """
    # The input crease pattern and angles
    data = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']
    
    # 1. Parse the input data
    angles = data[0::2]
    creases_in = data[1::2]
    
    print("Analyzing the given crease pattern:")
    print(f"Angles: {angles}")
    print(f"Crease assignments: {creases_in}\n")

    # For a vertex to be flat-foldable, it must satisfy both Kawasaki's and 
    # Maekawa's theorems. We will check these conditions to find the number 
    # of valid assignments.

    # 2. Check Kawasaki's Theorem
    # This theorem depends only on the angles. The sum of alternating angles 
    # around the vertex must be equal (and therefore 180 degrees).
    alt_sum1 = sum(angles[0::2])
    alt_sum2 = sum(angles[1::2])
    
    print("Checking Kawasaki's Theorem (sum of alternating angles must be equal):")
    
    # Output the numbers and result for the first set of alternating angles
    equation1_str = ' + '.join(map(str, angles[0::2]))
    print(f"Sum of odd-indexed angles: {equation1_str} = {alt_sum1}")
    
    # Output the numbers and result for the second set of alternating angles
    equation2_str = ' + '.join(map(str, angles[1::2]))
    print(f"Sum of even-indexed angles: {equation2_str} = {alt_sum2}")

    # For flat-foldability, the two sums must be equal.
    if alt_sum1 != alt_sum2:
        print("\nResult: Kawasaki's Theorem is NOT satisfied as the sums are unequal.")
        print("Therefore, no assignment of creases can make this pattern flat-foldable.")
        print("\nThe total number of different assignments is 0.")
        return

    # If Kawasaki's theorem were satisfied, we would proceed to check Maekawa's.
    # 3. Check Maekawa's Theorem for all possible assignments.
    # Maekawa's Theorem: |#Mountain - #Valley| = 2.
    
    initial_M = creases_in.count('M')
    initial_V = creases_in.count('V')
    q_indices = [i for i, c in enumerate(creases_in) if c == '?']
    num_q = len(q_indices)
    
    valid_assignment_count = 0
    
    print("\nResult: Kawasaki's Theorem is satisfied. Now checking Maekawa's Theorem for all assignments.")
    
    # Iterate through all 2^num_q ways to assign 'M' or 'V' to '?'
    for p in itertools.product(['M', 'V'], repeat=num_q):
        total_M = initial_M + p.count('M')
        total_V = initial_V + p.count('V')
        
        if abs(total_M - total_V) == 2:
            valid_assignment_count += 1
            
    print(f"\nFound {valid_assignment_count} assignments that satisfy Maekawa's Theorem.")
    print(f"The total number of different assignments is {valid_assignment_count}.")

# Execute the solution
if __name__ == "__main__":
    solve_flat_folding()
