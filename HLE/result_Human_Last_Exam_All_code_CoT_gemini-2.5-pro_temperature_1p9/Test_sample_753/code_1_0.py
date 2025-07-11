import math

def solve():
    """
    Solves the problem by calculating the cardinality of the set Theta^-1(lambda)
    and printing its first 40 digits.
    """
    
    # Step 1: Determine n for m=3
    m = 3
    n = sum(k * (m + 1 - k) for k in range(1, m + 1))

    # Step 2: Determine the partition lambda and its cycle structure
    # For m=3, the partition lambda is (3^1, 2^2, 1^3)
    # This corresponds to a permutation in Sigma_n with:
    # c_3 = 1 (one cycle of length 3)
    # c_2 = 2 (two cycles of length 2)
    # c_1 = 3 (three cycles of length 1)
    cycle_counts = {1: 3, 2: 2, 3: 1}

    # Step 3: Calculate the size of the centralizer |Z(pi)|
    centralizer_size = 1
    centralizer_formula_str = []
    centralizer_calc_str = []
    
    # Iterate through k in sorted order for display
    for k in sorted(cycle_counts.keys()):
        c_k = cycle_counts[k]
        term = (k**c_k) * math.factorial(c_k)
        centralizer_size *= term
        centralizer_formula_str.append(f"({k}^{c_k} * {c_k}!)")
        centralizer_calc_str.append(f"({k**c_k} * {math.factorial(c_k)})")

    full_centralizer_formula = " * ".join(centralizer_formula_str)
    full_centralizer_calc = " * ".join(centralizer_calc_str)
    
    # Step 4: Calculate the cardinality
    n_factorial = math.factorial(n)
    numerator = n_factorial**2
    cardinality = numerator // centralizer_size

    # Step 5: Format and print the output
    print(f"For m = {m}, n is calculated as:")
    print(f"n = sum_{{k=1}}^{{{m}}} k*({m}+1-k) = {n}")
    
    print("\nThe partition lambda corresponds to a permutation with cycle structure:")
    print(f"c_1={cycle_counts[1]}, c_2={cycle_counts[2]}, c_3={cycle_counts[3]}")

    print("\nThe cardinality of the set is calculated using the formula (n!)^2 / |Z(pi)|.")
    
    print("\n1. Calculate the size of the centralizer |Z(pi)|:")
    print(f"|Z(pi)| = {full_centralizer_formula}")
    print(f"|Z(pi)| = {full_centralizer_calc}")
    print(f"|Z(pi)| = {centralizer_size}")

    print("\n2. Calculate the numerator (n!)^2:")
    print(f"(n!)^2 = ({n}!)**2 = {n_factorial}^2 = {numerator}")

    print("\n3. Compute the final cardinality:")
    print(f"|Theta^-1(lambda)| = {numerator} / {centralizer_size} = {cardinality}")

    print("\nThe first 40 digits of the cardinality are:")
    # Pad the result with leading zeros to make it 40 digits long
    final_digits = str(cardinality).zfill(40)
    print(final_digits)
    print(f"\n<<<{final_digits}>>>")

solve()