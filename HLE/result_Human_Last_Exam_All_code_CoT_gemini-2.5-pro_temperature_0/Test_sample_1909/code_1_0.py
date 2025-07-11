def solve_ratio_limit(n_max):
    """
    Calculates the expected number of remaining numbers, E_n,
    and the ratio E_n/n up to n_max.
    """
    if n_max < 2:
        print("n_max must be at least 2.")
        return

    # Initialize E with E_0 and E_1
    e = [0.0, 1.0]

    # Use the simplified recurrence E_n = ((n-2)E_{n-1} + 2E_{n-2}) / (n-1)
    # which is valid for n >= 2.
    # For n=2: E_2 = (0*E_1 + 2*E_0)/1 = 0
    for n in range(2, n_max + 1):
        e_n = ((n - 2) * e[n - 1] + 2 * e[n - 2]) / (n - 1)
        e.append(e_n)

    ratio = e[n_max] / n_max
    
    # The problem asks for the equation for the final step of the recurrence
    # E_n = ((n-2)*E_{n-1} + 2*E_{n-2}) / (n-1)
    # We will print the calculation for E_n_max
    
    n = n_max
    e_n_minus_1 = e[n-1]
    e_n_minus_2 = e[n-2]
    e_n = e[n]
    
    print(f"Calculation for E_{n}:")
    print(f"E_{n} = (({n}-2) * E_{{{n}-1}} + 2 * E_{{{n}-2}}) / ({n}-1)")
    print(f"E_{n} = (({n-2}) * {e_n_minus_1:.6f} + 2 * {e_n_minus_2:.6f}) / {n-1}")
    print(f"E_{n} = {e_n:.6f}")
    print("-" * 20)
    print(f"The ratio E_{n}/n for n={n} is:")
    print(f"{e_n:.6f} / {n} = {ratio:.6f}")
    print("-" * 20)
    print(f"The limit of the ratio as n approaches infinity is e^(-2).")
    print(f"e^(-2) is approximately {0.135335}")


# Set a large value for n to approximate the limit
n_limit = 1000
solve_ratio_limit(n_limit)
