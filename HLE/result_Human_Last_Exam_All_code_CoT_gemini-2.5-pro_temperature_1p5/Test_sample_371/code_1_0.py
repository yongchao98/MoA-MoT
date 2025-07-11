def solve_and_print():
    """
    This function implements the optimal solution to find the k-th element of S_n,
    and prints the calculation steps as requested.
    """
    # Per the prompt, we use n and k values larger than a typical INT_MAX/2.
    # On a 64-bit system, sys.maxsize is 2**63 - 1. We'll use k slightly larger
    # than half of that. Note that the value of n is irrelevant in the optimal solution.
    n = 2**63
    k = (2**63 - 1) // 2 + 12345

    # The optimal solution relies on the insight that S_n[k] = v_2(k+1) + 1,
    # where v_2(x) is the 0-indexed position of the least significant bit (LSB) of x.
    # This is calculated efficiently as ((k+1) & -(k+1)).bit_length().
    x = k + 1

    # Isolate the value of the LSB, which is always a power of 2.
    lsb_val = x & -x

    # For any power of two, y, its bit_length() gives the 1-based index of the set bit.
    result = lsb_val.bit_length()

    # The prompt asks to output each number in the final equation.
    # We will show the formula and substitute the values step-by-step.
    print(f"To find S_n[k] for n={n}, k={k}:")
    print(f"The value is given by the formula: S_n[k] = bit_length( (k+1) & -(k+1) )")
    print(f"Substituting k: S_n[{k}] = bit_length( ({k}+1) & -({k}+1) )")
    print(f"Calculating k+1: S_n[{k}] = bit_length( {x} & -{x} )")
    print(f"Calculating the LSB: S_n[{k}] = bit_length( {lsb_val} )")
    print(f"Final result: S_n[{k}] = {result}")

solve_and_print()