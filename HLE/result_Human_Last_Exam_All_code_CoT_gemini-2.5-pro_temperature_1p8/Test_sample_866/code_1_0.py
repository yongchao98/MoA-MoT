import itertools

def solve_flat_folding():
    """
    Solves the flat-folding problem by checking the necessary theorems
    and counting the number of valid crease assignments.
    """
    # The input describes a single vertex with 5 angles and 5 creases.
    # Pattern: [60,M,30,?,50,?,70,V,150,?]
    angles = [60, 30, 50, 70, 150]
    creases_template = ['M', '?', '?', 'V', '?']

    print("Analyzing the crease pattern for flat-foldability.")
    print("--------------------------------------------------")

    # --- Step 1: Check Kawasaki's Theorem (Geometric Condition) ---
    print("Step 1: Checking Kawasaki's Theorem.")
    print("This theorem states that the sums of alternate angles must be equal.")
    sum_odd_angles = sum(angles[0::2])
    sum_even_angles = sum(angles[1::2])
    
    print(f"The angles are: {angles}")
    print(f"Sum of odd-indexed angles: { ' + '.join(map(str, angles[0::2])) } = {sum_odd_angles}")
    print(f"Sum of even-indexed angles: { ' + '.join(map(str, angles[1::2])) } = {sum_even_angles}")

    if sum_odd_angles != sum_even_angles:
        print("\nResult: Kawasaki's Theorem is NOT satisfied.")
        print("This is a fundamental geometric requirement for the paper to lie flat.")
        print("As this condition fails, no M/V assignment can make the pattern flat-foldable.")
        print("--------------------------------------------------")
        print("Final Equation: 0 = 0")
        print("\nTotal valid assignments: 0")
        return

    # --- Step 2: Check Maekawa's Theorem (Crease Assignment Condition) ---
    print("\nStep 2: Checking Maekawa's Theorem.")
    print("This theorem states |#Mountain - #Valley| = 2.")
    num_creases = len(creases_template)
    print(f"The pattern has {num_creases} creases in total.")
    
    if num_creases % 2 != 0:
        print("Note: The total number of creases is odd. Maekawa's Theorem can only be satisfied if the total number of creases is even.")
        print("We expect to find no valid assignments.")

    # Find unknown creases and iterate through all possibilities
    unknown_indices = [i for i, c in enumerate(creases_template) if c == '?']
    num_unknown = len(unknown_indices)
    base_m_count = creases_template.count('M')
    base_v_count = creases_template.count('V')
    
    all_assignments = list(itertools.product(['M', 'V'], repeat=num_unknown))
    validity_counts = []

    print(f"\nThere are {num_unknown} unassigned creases, giving {len(all_assignments)} possible combinations to test:")
    for p in all_assignments:
        m_count = base_m_count + p.count('M')
        v_count = base_v_count + p.count('V')
        
        # Determine if this assignment is valid
        is_valid = 1 if abs(m_count - v_count) == 2 else 0
        validity_counts.append(is_valid)
    
    # --- Step 3: Calculate Final Result and Display Equation ---
    print("\n--------------------------------------------------")
    print("Step 3: Final Calculation.")
    print("The final count is the sum of outcomes (1 for a valid assignment, 0 for invalid).")

    total_valid = sum(validity_counts)
    equation_str = " + ".join(map(str, validity_counts))
    
    print(f"Final Equation: {equation_str} = {total_valid}")
    print(f"\nTotal number of different assignments which make the pattern flat foldable is: {total_valid}")


# Run the solver
solve_flat_folding()
<<<0>>>