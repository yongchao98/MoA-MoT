import math

def solve_augustus_records():
    """
    Solves the problem by establishing and solving a system of equations
    based on the provided constraints.
    """
    
    # Step 1 & 2: Establish knowns and interpret constraints
    # 'A' (dual-named) + 'C' (full title) make up the 1/4 of readable records with "Augustus" and "Caesar".
    # Total readable records = 720 * (5/6) = 600.
    # So, A + C = 600 / 4 = 150.
    A_plus_C = 150

    # Constraint for 'B' (single-variant) is sqrt(B) = (lost records)/4.
    # Assuming "lost records" is a riddle for the 80 "Caesar" records to avoid contradiction.
    caesar_only_records = 80
    B_sqrt = caesar_only_records / 4
    B = int(B_sqrt ** 2)

    # Step 3: Solve the system of equations by iterating
    # We have A + C = 150 and sqrt(A) + 3*L = C.
    # This implies 3*L = 150 - A - sqrt(A).
    # We need to find an integer L, and A must be a perfect square.
    
    possible_solutions = []
    # A must be a perfect square less than 150
    for y in range(1, int(math.sqrt(A_plus_C)) + 1):
        A_candidate = y * y
        # Check if 150 - A - sqrt(A) is a positive multiple of 3
        if (A_plus_C - A_candidate - y) >= 0 and (A_plus_C - A_candidate - y) % 3 == 0:
            L_candidate = (A_plus_C - A_candidate - y) // 3
            possible_solutions.append({'L': L_candidate, 'A': A_candidate})

    # Step 4: Use the divisibility rule to find the unique solution.
    # Sum of found numbers (L+A+B+C) must be divisible by the number of naming patterns.
    # Sum = L + A + B + C = L + (A+C) + B = L + 150 + 400 = L + 550.
    # Denominator: distinct ways Augustus is named:
    # 1. "Augustus"+"Caesar" (A), 2. "Octavianus" (B), 3. "Augustus" (B), 4. "Imp. Cae. Aug." (C)
    # Divisor = 4.
    divisor = 4
    
    final_solution = None
    for sol in possible_solutions:
        L_val = sol['L']
        sum_of_numbers = L_val + A_plus_C + B
        if sum_of_numbers % divisor == 0:
            final_solution = sol
            break
            
    # Step 5: Calculate final values based on the unique solution
    L = final_solution['L']
    A = final_solution['A']
    C = A_plus_C - A
    
    print("--- Calculated Values ---")
    print(f"Lost records (L): {L}")
    print(f"Dual-named documents (A): {A}")
    print(f"Single-variant documents (B): {B}")
    print(f"Full title documents (C): {C}")
    print(f"Number of distinct naming patterns for the final division: {divisor}")
    print("\n--- Final Calculation ---")
    
    product = L * A * B * C
    
    print(f"Equation: ({L} * {A} * {B} * {C}) / {divisor}")
    
    answer = product / divisor
    final_answer = math.floor(answer / 1000)
    
    print(f"Result of division: {answer}")
    print(f"Floor of result / 1000: {final_answer}")
    
    return final_answer

final_answer_val = solve_augustus_records()
print(f"<<<{final_answer_val}>>>")