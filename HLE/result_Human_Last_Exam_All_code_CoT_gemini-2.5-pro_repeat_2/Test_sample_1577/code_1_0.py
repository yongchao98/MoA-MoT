import math

def calculate_toric_code_degeneracy(n, m):
    """
    Calculates the ground space degeneracy of the toric code with n smooth and m rough holes.

    Args:
        n (int): The number of smooth holes.
        m (int): The number of rough holes.
    """
    if n < 0 or m < 0 or not isinstance(n, int) or not isinstance(m, int):
        print("Error: Number of holes must be a non-negative integer.")
        return

    # The formula for the number of logical qubits k is k = n + m - 2.
    # This is based on a model counting independent loop operators and is most
    # applicable when n>=1 and m>=1.
    # The ground space degeneracy (GSD) is 2^k.
    
    print(f"For a toric code with n = {n} smooth holes and m = {m} rough holes:")
    
    if n + m < 2:
        # The formula n+m-2 is not applicable for cases with 0 or 1 total holes.
        # GSD for 0 holes (sphere) is 1.
        # GSD for 1 hole (disk) is 1.
        if n + m == 0:
            exponent_str = "0 (by definition for a sphere)"
            degeneracy = 1
            print(f"The ground space degeneracy is {degeneracy}.")
        elif n + m == 1:
            exponent_str = "0 (by definition for a disk)"
            degeneracy = 1
            print(f"The ground space degeneracy is {degeneracy}.")
        else: # Should not be reached
             print("Invalid case.")
    else:
        exponent = n + m - 2
        # Use math.pow to handle potentially large numbers, will return float
        degeneracy = int(math.pow(2, exponent))
        print(f"The degeneracy is given by the formula: 2^(n + m - 2)")
        print(f"GSD = 2^({n} + {m} - 2) = 2^{exponent} = {degeneracy}")


# Example values for n and m. You can change these to see other results.
n_smooth_holes = 3
m_rough_holes = 4

calculate_toric_code_degeneracy(n_smooth_holes, m_rough_holes)