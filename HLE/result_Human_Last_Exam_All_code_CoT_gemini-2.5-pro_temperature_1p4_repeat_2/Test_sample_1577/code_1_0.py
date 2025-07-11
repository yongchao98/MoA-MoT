def calculate_toric_code_gsd(n, m):
    """
    Calculates the ground space degeneracy (GSD) of the toric code with n smooth holes and m rough holes.

    Args:
        n (int): The number of smooth holes (must be >= 1).
        m (int): The number of rough holes (must be >= 1).
    """
    if n < 1 or m < 1:
        print("This formula 2^(n+m-2) is valid for n >= 1 and m >= 1.")
        # The more general formula gives 2**(max(0, n-1) + max(0, m-1))
        # For this problem, we stick to the formula from the likely answer choice.
        return

    # The number of logical qubits is k = (n-1) + (m-1) = n+m-2
    k = n + m - 2
    
    # The GSD is 2^k
    gsd = 2**k

    print("The ground space degeneracy (GSD) is given by the formula: 2^(n + m - 2)")
    print(f"For n = {n} and m = {m}:")
    # Using f-string to explicitly show the calculation steps in the final equation.
    print(f"GSD = 2^({n} + {m} - 2) = 2^{k} = {gsd}")

# Example usage with n=2 smooth holes and m=3 rough holes.
calculate_toric_code_gsd(2, 3)
