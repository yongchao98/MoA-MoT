def solve_endomorphisms():
    """
    Calculates the number of equivalence classes of endomorphisms on a set of size 4.

    This is equivalent to finding the number of non-isomorphic functional graphs
    on 4 vertices, which we denote as a_4.
    """
    n = 4
    
    # Step 1: Pre-computed values for T_n, the number of connected non-isomorphic
    # functional graphs on n vertices. (OEIS A000273)
    T = {1: 1, 2: 2, 3: 4, 4: 9}
    
    print("This program calculates the number of equivalence classes of endomorphisms on a set of size 4.")
    print("Let a_n be the number of such classes for a set of size n.")
    print("Let T_n be the number of *connected* such graphs.\n")
    print("Known values for T_n:")
    for i in range(1, n + 1):
        print(f"T_{i} = {T[i]}")
    print("-" * 30)

    # Step 2: Calculate intermediate coefficients b_i, where b_i = sum_{d|i} d * T_d.
    b = {}
    print("Step 1: Calculate intermediate coefficients b_i")
    for i in range(1, n + 1):
        b[i] = 0
        divisors = [d for d in range(1, i + 1) if i % d == 0]
        sum_str_list = []
        for d in divisors:
            b[i] += d * T[d]
            sum_str_list.append(f"{d}*T[{d}]")
        print(f"b_{i} = {' + '.join(sum_str_list)} = {b[i]}")
    print("-" * 30)

    # Step 3: Use the recurrence relation to find a_n.
    # Recurrence: n * a_n = sum_{i=1 to n} a_{n-i} * b_i
    # Initialize a_0 = 1 (for the empty graph on 0 vertices).
    a = {0: 1}
    print("Step 2: Calculate a_n using the recurrence relation")
    print("a_0 = 1")
    for i in range(1, n + 1):
        numerator = 0
        # Build the numerator of the recurrence relation
        for j in range(1, i + 1):
            numerator += a[i - j] * b[j]
        
        # Store the result
        a[i] = numerator // i
        
        # Format the output to show the calculation step-by-step
        num_symbol_str = " + ".join([f"a_{i-j}*b_{j}" for j in range(1, i + 1)])
        num_value_str = " + ".join([f"{a[i-j]}*{b[j]}" for j in range(1, i + 1)])
        
        print(f"\nCalculating a_{i}:")
        print(f"{i} * a_{i} = {num_symbol_str}")
        print(f"{i} * a_{i} = {num_value_str}")
        print(f"{i} * a_{i} = {numerator}")
        print(f"a_{i} = {numerator} / {i} = {a[i]}")

    print("-" * 30)
    print(f"The number of elements of E represented by F is a_4 = {a[4]}.")

solve_endomorphisms()