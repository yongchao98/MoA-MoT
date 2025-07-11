import math

def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent alpha (α) using the first-order
    epsilon expansion in renormalization group theory.
    """
    # 1. Define physical parameters
    # The spatial dimension of the system.
    d = 3
    # The upper critical dimension, above which mean-field theory is valid. For O(n) models, it is 4.
    d_c = 4
    # The number of components of the order parameter. We assume the Ising universality class, where n=1.
    n = 1

    # 2. Calculate the expansion parameter, epsilon
    epsilon = d_c - d

    # 3. Use the first-order epsilon expansion formula to calculate alpha
    # Calculate the numerator and denominator for the fractional result
    numerator = (4 - n) * epsilon
    denominator = 2 * (n + 8)

    # Calculate the final decimal value
    alpha_value = numerator / denominator
    
    # --- Output Section ---
    print("This script calculates the specific heat scaling exponent α using the ε-expansion.")
    print("-" * 60)
    print("1. Parameters for the calculation:")
    print(f"   - Spatial dimension (d): {d}")
    print(f"   - Upper critical dimension (d_c): {d_c}")
    print(f"   - Order parameter components (n, for Ising model): {n}")
    print("-" * 60)

    print("2. The expansion parameter ε is calculated as:")
    print(f"   ε = d_c - d = {d_c} - {d} = {epsilon}")
    print("-" * 60)
    
    print("3. The formula for α to first order in ε is:")
    print("   α = (4 - n) / (2 * (n + 8)) * ε")
    print("-" * 60)

    print("4. Substituting the values into the equation:")
    # The user requested to see each number in the final equation.
    print(f"   α = (4 - {n}) / (2 * ({n} + 8)) * {epsilon}")
    print(f"   α = ({4 - n}) / (2 * {n + 8})")
    print(f"   α = {numerator} / {denominator}")
    print("-" * 60)
    
    # Simplify the fraction for a more elegant representation (e.g., 3/18 -> 1/6)
    common_divisor = math.gcd(numerator, denominator)
    simplified_num = numerator // common_divisor
    simplified_den = denominator // common_divisor
    
    print("5. Final Result:")
    print(f"   The scaling exponent α is {simplified_num}/{simplified_den}, which is approximately {alpha_value}.")

if __name__ == "__main__":
    calculate_alpha_exponent()
    # The final numerical answer in the requested format
    final_alpha = 1/6
    # print(f"\n<<<{final_alpha}>>>") # This is a comment for the final answer extraction