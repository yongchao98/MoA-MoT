import itertools

def solve_crease_pattern():
    """
    Analyzes a given partial crease pattern to find the number of valid
    flat-foldable assignments.
    """
    raw_input = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']
    
    # 1. Parse input into angles and creases
    angles = [item for item in raw_input if isinstance(item, int)]
    creases = [item for item in raw_input if isinstance(item, str)]
    
    print(f"Input Crease Pattern: {raw_input}\n")
    print("For a vertex to be flat-foldable, it must satisfy several conditions.")
    
    # 2. Check Kawasaki's Theorem
    # For a fixed sequence of angles, the sum of alternating angles must be 180.
    sum_odd = sum(angles[i] for i in range(0, len(angles), 2))
    sum_even = sum(angles[i] for i in range(1, len(angles), 2))
    
    print("--- Condition 1: Kawasaki's Theorem (Angles) ---")
    print(f"The sequence of angles is: {angles}")
    print(f"The sum of angles must be 360 degrees. {sum(angles)} = 360. This is satisfied.")
    print("The alternating sums of angles must be equal (and thus 180).")
    print(f"Sum of odd-positioned angles: {angles[0]} + {angles[2]} + {angles[4]} = {sum_odd}")
    print(f"Sum of even-positioned angles: {angles[1]} + {angles[3]} = {sum_even}")

    kawasaki_ok = (sum_odd == 180 and sum_even == 180)
    if kawasaki_ok:
        print("Result: Kawasaki's Theorem is satisfied.\n")
    else:
        print(f"Result: Since {sum_odd} != {sum_even}, Kawasaki's Theorem is NOT satisfied for this angle sequence.\n")

    # 3. Check Maekawa's Theorem for all assignments
    # The number of mountain (M) and valley (V) creases must differ by 2. |M - V| = 2.
    print("--- Condition 2: Maekawa's Theorem (Creases) ---")
    print("We will now check all possible assignments for the '?' creases.")
    
    q_indices = [i for i, c in enumerate(creases) if c == '?']
    num_q = len(q_indices)
    valid_assignments_count = 0
    
    # Iterate through all possible assignments for '?'
    assignment_options = list(itertools.product(['M', 'V'], repeat=num_q))

    for i, p in enumerate(assignment_options):
        temp_creases = list(creases)
        # Apply the current assignment combination
        for j in range(num_q):
            temp_creases[q_indices[j]] = p[j]
        
        m_count = temp_creases.count('M')
        v_count = temp_creases.count('V')
        
        print(f"\nChecking Assignment #{i+1} of {len(assignment_options)}:")
        print(f"Crease assignment: {temp_creases}")
        print(f"Equation: |Number of Mountains - Number of Valleys| = 2")
        print(f"Calculation: |{m_count} - {v_count}| = {abs(m_count - v_count)}")
        
        maekawa_ok = abs(m_count - v_count) == 2
        if maekawa_ok:
            print("Result: Maekawa's condition is SATISFIED.")
            # A fully valid assignment must satisfy both theorems
            if kawasaki_ok:
                valid_assignments_count += 1
        else:
            print("Result: Maekawa's condition is NOT satisfied.")
    
    print("\n--- Final Count ---")
    print("An assignment is valid only if BOTH Kawasaki's and Maekawa's theorems are satisfied.")
    print(f"Total number of assignments satisfying Maekawa's Theorem: {valid_assignments_count}")
    print(f"Is Kawasaki's Theorem satisfied by the angles? {'Yes' if kawasaki_ok else 'No'}")
    
    final_answer = 0
    if kawasaki_ok:
        final_answer = valid_assignments_count

    print(f"\nTherefore, the total number of different assignments which will make the pattern flat fordable is {final_answer}.")

solve_crease_pattern()
print("<<<0>>>")