import math

# --- Problem Parameters ---
# Size of the grid is n x m
n = 14
m = 14

# The rule for infection requires at least k infected neighbors
k = 3

# --- Calculation based on the formula from Benevides & Przykucki (2013) ---

# First, we calculate the intermediate values n' and m'
k_minus_1 = k - 1
n_minus_1 = n - 1
m_minus_1 = m - 1

n_prime = math.floor(n_minus_1 / k_minus_1)
m_prime = math.floor(m_minus_1 / k_minus_1)

# Now, we apply the main formula: m_k(n, m) = (k-1) * (n' + m') + 2
minimum_sites = k_minus_1 * (n_prime + m_prime) + 2

# --- Output the reasoning and the result ---
print("This problem is a case of bootstrap percolation on a grid.")
print(f"We are considering an {n}x{m} grid with a {k}-neighbor rule.")
print("\nAccording to the formula from Benevides & Przykucki (2013), the minimum number of initial sites is calculated as:")
print("m_k(n, m) = (k-1) * (n' + m') + 2")
print("where n' = floor((n-1)/(k-1)) and m' = floor((m-1)/(k-1)).")

print("\nStep 1: Calculate n' and m'")
print(f"k-1 = {k} - 1 = {k_minus_1}")
print(f"n' = floor(({n}-1)/{k_minus_1}) = floor({n_minus_1}/{k_minus_1}) = {n_prime}")
print(f"m' = floor(({m}-1)/{k_minus_1}) = floor({m_minus_1}/{k_minus_1}) = {m_prime}")

print("\nStep 2: Substitute these values into the main formula")
print(f"Result = {k_minus_1} * ({n_prime} + {m_prime}) + 2")
print(f"Result = {k_minus_1} * ({n_prime + m_prime}) + 2")
print(f"Result = {k_minus_1 * (n_prime + m_prime)} + 2")
print(f"Result = {minimum_sites}")

print("\nTherefore, the minimum number of initially infected sites is:")
print(minimum_sites)