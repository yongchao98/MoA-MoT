import math

def combinations(n, k):
    """
    Computes the binomial coefficient C(n, k) = n! / (k! * (n-k)!).
    """
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    # Use symmetry C(n, k) = C(n, n-k)
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        # Perform multiplication before division to maintain precision
        res = res * (n - i) // (i + 1)
    return res

def laguerre_L_at_1(n):
    """
    Calculates the value of the generalized Laguerre polynomial
    L_{n+1}^{(3n-1)}(1) using its series definition.
    k = n + 1, alpha = 3n - 1
    """
    k = n + 1
    alpha = 3 * n - 1
    k_plus_alpha = k + alpha
    
    total_sum = 0
    for i in range(k + 1):
        # L_k^(alpha)(x) = sum_{i=0 to k} C(k+alpha, k-i) * (-x)^i / i!
        # Here x=1
        term = combinations(k_plus_alpha, k - i) * ((-1)**i) / math.factorial(i)
        total_sum += term
    return total_sum

def Mz1(n):
    """
    Calculates the magnetization M_z(1) for a given number of spins, n.
    """
    # The number of spins n must be a positive integer.
    if n <= 0 or not isinstance(n, int):
        return float('inf')

    # Calculate the value of the Laguerre polynomial term
    laguerre_val = laguerre_L_at_1(n)
    
    # Calculate the prefactor term
    numerator = ((-1)**n) * (n + 1) * (n**(-n))
    denominator = (math.pi / 2)**n
    prefactor = numerator / denominator
    
    # Calculate the final M_z(1) value
    result = prefactor * laguerre_val
    return result, prefactor, laguerre_val

# --- Main execution ---

# Search for the n that minimizes M_z(1) in a practical range for n.
min_Mz = float('inf')
n_min = -1
best_prefactor = 0
best_laguerre = 0

# A search range up to 15 should be sufficient to find the minimum.
for n_val in range(1, 16):
    current_Mz, _, _ = Mz1(n_val)
    if current_Mz < min_Mz:
        min_Mz, best_prefactor, best_laguerre = Mz1(n_val)
        n_min = n_val

# --- Output the final answer ---
print(f"The minimum magnetization M_z(1) is found for n = {n_min}.")
print("\nThe final value is calculated as follows:")
print(f"n = {n_min}")
print(f"Prefactor = ((-1)^{n_min} * ({n_min} + 1) * {n_min}^-{n_min}) / (pi / 2)^{n_min}")
print(f"Prefactor value = {best_prefactor}")
print(f"Laguerre polynomial term L_{{{n_min}+1}}^{{(3*{n_min}-1)}}(1) value = {best_laguerre}")
print(f"\nM_z(1) = Prefactor * Laguerre_term")
print(f"M_z(1) = {best_prefactor} * {best_laguerre}")
print(f"Final minimum magnetization M_z(1) = {min_Mz}")

<<< -11.69837330752528 >>>