def solve_lie_algebra_problem():
    """
    Calculates the number of non-compact positive roots for a real form of C8.

    The real form is specified by a Vogan diagram, which tells us which simple
    roots are non-compact. For a positive root expressed as a sum of simple
    roots, it is non-compact if the sum of coefficients of non-compact simple
    roots is odd.

    The script iterates through all positive roots of C8, determines their
    coefficients in the simple root basis, and checks the non-compactness
    condition.
    """
    N = 8
    # The set of indices of non-compact simple roots, from the Vogan diagram
    # W-B-W-B-B-W-B==B corresponds to alpha_1, alpha_3, alpha_6 being non-compact.
    non_compact_simple_indices = {1, 3, 6}

    counts = {
        "type1": 0, # e_i - e_j
        "type2": 0, # e_i + e_j, j < N
        "type3": 0, # e_i + e_N
        "type4": 0, # 2*e_i, i < N
        "type5": 0, # 2*e_N
    }

    # Type 1: roots of the form e_i - e_j for 1 <= i < j <= N
    # Expansion: sum_{k=i}^{j-1} alpha_k
    for i in range(1, N + 1):
        for j in range(i + 1, N + 1):
            coeffs = [0] * N
            for k in range(i, j):
                coeffs[k - 1] = 1
            
            parity_sum = sum(coeffs[k - 1] for k in non_compact_simple_indices)
            if parity_sum % 2 != 0:
                counts["type1"] += 1

    # Type 2: roots of the form e_i + e_j for 1 <= i < j < N
    # Expansion: sum_{k=i}^{j-1} alpha_k + 2*sum_{k=j}^{N-1} alpha_k + alpha_N
    for i in range(1, N):
        for j in range(i + 1, N):
            coeffs = [0] * N
            for k in range(i, j):
                coeffs[k - 1] = 1
            for k in range(j, N):
                coeffs[k - 1] = 2
            coeffs[N - 1] = 1
            
            parity_sum = sum(coeffs[k - 1] for k in non_compact_simple_indices)
            if parity_sum % 2 != 0:
                counts["type2"] += 1

    # Type 3: roots of the form e_i + e_N for 1 <= i < N
    # Expansion: sum_{k=i}^{N-1} alpha_k + alpha_N
    for i in range(1, N):
        coeffs = [0] * N
        for k in range(i, N):
            coeffs[k - 1] = 1
        coeffs[N - 1] = 1 # This is redundant from the loop, but for clarity
        
        parity_sum = sum(coeffs[k - 1] for k in non_compact_simple_indices)
        if parity_sum % 2 != 0:
            counts["type3"] += 1

    # Type 4: roots of the form 2*e_i for 1 <= i < N
    # Expansion: 2*sum_{k=i}^{N-1} alpha_k + alpha_N
    for i in range(1, N):
        coeffs = [0] * N
        for k in range(i, N):
            coeffs[k - 1] = 2
        coeffs[N - 1] = 1
        
        parity_sum = sum(coeffs[k - 1] for k in non_compact_simple_indices)
        if parity_sum % 2 != 0:
            counts["type4"] += 1

    # Type 5: root of the form 2*e_N
    # Expansion: alpha_N
    coeffs = [0] * N
    coeffs[N - 1] = 1
    parity_sum = sum(coeffs[k - 1] for k in non_compact_simple_indices)
    if parity_sum % 2 != 0:
        counts["type5"] += 1
        
    total_non_compact = sum(counts.values())

    print(f"Number of non-compact positive roots of type e_i - e_j: {counts['type1']}")
    print(f"Number of non-compact positive roots of type e_i + e_j (j < 8): {counts['type2']}")
    print(f"Number of non-compact positive roots of type e_i + e_8: {counts['type3']}")
    print(f"Number of non-compact positive roots of type 2*e_i (i < 8): {counts['type4']}")
    print(f"Number of non-compact positive roots of type 2*e_8: {counts['type5']}")
    
    equation_parts = [str(c) for c in counts.values()]
    equation_str = " + ".join(equation_parts)
    print(f"Total number of non-compact positive roots = {equation_str} = {total_non_compact}")

solve_lie_algebra_problem()