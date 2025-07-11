import math

def solve_crease_pattern():
    """
    Calculates the number of valid flat-foldable assignments for a given
    crease pattern at a single vertex.
    """
    # The input crease pattern provided by the user.
    pattern_str = "[60,M,30,?,50,?,70,V,150,?]"
    
    # Parse the input string into lists of angles and crease specifications.
    # Note: Using eval() is generally unsafe with untrusted input, but acceptable here
    # as the input is fixed in the problem description.
    try:
        pattern = eval(pattern_str)
        angles = pattern[::2]
        creases_spec = pattern[1::2]
    except (SyntaxError, TypeError):
        print("Error: Invalid input format.")
        return

    n = len(angles)
    if n != len(creases_spec):
        print("Error: The pattern must have an equal number of angles and creases.")
        return

    # A basic check: for a vertex on a flat plane, angles must sum to 360.
    if not math.isclose(sum(angles), 360):
        print(f"Sum of angles is {sum(angles)}, which is not 360. Cannot form a flat vertex.")
        print(0)
        return

    # Condition 1: Check Kawasaki's Theorem.
    # The alternating sum of angles must be zero for the vertex to be flat-foldable.
    alternating_sum = sum(angle * ((-1)**i) for i, angle in enumerate(angles))
    
    print("Step 1: Checking Kawasaki's Theorem (alternating sum of angles must be 0)")
    equation_parts = []
    for i, angle in enumerate(angles):
        if i > 0:
            equation_parts.append("+" if i % 2 == 0 else "-")
        equation_parts.append(str(angle))
    
    equation_str = " ".join(equation_parts)
    print(f"Equation: {equation_str} = {alternating_sum}")

    if not math.isclose(alternating_sum, 0):
        print("Result: Kawasaki's Theorem is NOT satisfied.")
        print("\nSince Kawasaki's theorem is a necessary condition, no assignment can make this pattern flat-foldable.")
        print("Total number of different assignments: 0")
        return

    # Condition 2: Check Maekawa's Theorem.
    # The number of mountain (N_M) and valley (N_V) folds must differ by 2.
    print("\nStep 2: Checking Maekawa's Theorem (|N_M - N_V| = 2)")
    print(f"Total number of creases, n = {n}")
    print(f"This means N_M + N_V = {n}.")
    
    # Check if integer solutions exist for N_M and N_V
    # N_M - N_V = 2  => 2*N_M = n+2
    # N_V - N_M = 2  => 2*N_V = n+2
    if (n + 2) % 2 != 0:
        print("Result: Maekawa's Theorem can NOT be satisfied.")
        print("Reason: For |N_M - N_V| = 2 and N_M + N_V = n, N_M and N_V cannot be integers if n is odd.")
        print("\nSince Maekawa's theorem is a necessary condition, no assignment can make this pattern flat-foldable.")
        print("Total number of different assignments: 0")
        return

    # If both preconditions were met, we would count combinations.
    # This part of the code will not be reached for the given input.
    print("\nBoth Kawasaki's and Maekawa's theorems are satisfiable in principle.")
    print("Step 3: Counting valid assignments for '?'")
    
    num_unknowns = creases_spec.count('?')
    known_M = creases_spec.count('M')
    known_V = creases_spec.count('V')
    
    total_valid_assignments = 0
    
    # Case A: N_M - N_V = 2
    target_M_A = (n + 2) // 2
    needed_M_A = target_M_A - known_M
    if 0 <= needed_M_A <= num_unknowns:
        # The rest of unknowns must be V
        needed_V_A = num_unknowns - needed_M_A
        if needed_V_A + known_V == n - target_M_A:
            total_valid_assignments += math.comb(num_unknowns, needed_M_A)
            
    # Case B: N_V - N_M = 2
    target_M_B = (n - 2) // 2
    needed_M_B = target_M_B - known_M
    if 0 <= needed_M_B <= num_unknowns:
       # The rest of unknowns must be V
        needed_V_B = num_unknowns - needed_M_B
        if needed_V_B + known_V == n - target_M_B:
            # Prevent double counting if n+2 and n-2 are the same modulo something
            if target_M_A != target_M_B:
                total_valid_assignments += math.comb(num_unknowns, needed_M_B)

    print(f"\nTotal number of different assignments: {total_valid_assignments}")


solve_crease_pattern()