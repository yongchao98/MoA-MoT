import math

def calculate_toric_code_gsd(n, m):
    """
    Calculates the Ground Space Degeneracy (GSD) of the toric code
    on a surface with n smooth and m rough holes.

    This calculation assumes the general case where n > 0 and m > 0.

    Args:
        n (int): The number of smooth holes (must be > 0).
        m (int): The number of rough holes (must be > 0).
    """
    if n <= 0 or m <= 0:
        print("This formula is for the general case where n > 0 and m > 0.")
        # We can use the more general formula for these cases if needed,
        # but based on the problem's multiple-choice options,
        # we focus on the case n>0, m>0.
        k = n + m - 2 + (1 if n == 0 else 0) + (1 if m == 0 else 0)
        if k < 0: # This should not happen for n,m >= 0
            k = 0
        gsd = int(math.pow(2, k))
        print(f"For the specific case n={n}, m={m}, the degeneracy is {gsd}.")
        return

    # In the general case (n>0, m>0), the number of logical qubits k is n + m - 2.
    k = n + m - 2
    
    # The Ground Space Degeneracy (GSD) is 2^k.
    gsd = int(math.pow(2, k))

    print(f"For a toric code with n = {n} smooth holes and m = {m} rough holes:")
    print("The number of logical qubits k is given by the formula:")
    print("k = n + m - 2")
    print(f"k = {n} + {m} - 2")
    print(f"k = {k}")
    print("\nThe Ground Space Degeneracy (GSD) is 2^k:")
    print(f"GSD = 2^{k}")
    print(f"GSD = {gsd}")


# Example usage with n=3 smooth holes and m=2 rough holes.
calculate_toric_code_gsd(n=3, m=2)