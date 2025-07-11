import math

def main():
    """
    This script calculates the total mass for a given q and n.
    """
    # --- Parameters ---
    # You can modify these values to solve for your specific case.
    # q must be a prime power greater than 1.
    # n must be an integer greater than or equal to 1.
    q = 2
    n = 2

    # --- Functions ---
    def dedekind_zeta_function(k, q_val):
        """
        Calculates the value of the Dedekind zeta function Z(k) for the
        rational function field F_q(t). The formula is Z(s) = 1/((1-q^{-s})(1-q^{1-s})).
        """
        if k == 1:
            return float('inf')
        
        term1 = 1 - q_val**(-k)
        term2 = 1 - q_val**(1 - k)
        return 1 / (term1 * term2)

    # --- Calculation and Output ---
    print(f"Computing the total mass for q = {q} and n = {n}.")
    print("Formula: Total Mass = (q / (q - 1)) * PRODUCT_{k=2 to n} Z(k)")
    print("where Z(k) = 1 / ( (1 - q**(-k)) * (1 - q**(1-k)) )")
    print("-" * 40)

    # Calculate the product of zeta function values
    if n == 1:
        zeta_product = 1.0
        print("For n = 1, the product term is empty, which evaluates to 1.")
    else:
        zeta_product = 1.0
        zeta_values_str = []
        print("Step 1: Calculate Z(k) for k from 2 to n.")
        for k in range(2, n + 1):
            z_k_val = dedekind_zeta_function(k, q)
            print(f"Z({k}) = 1 / ((1 - {q}**(-{k}))*(1 - {q}**(1-{k}))) = {z_k_val}")
            zeta_values_str.append(str(z_k_val))
            zeta_product *= z_k_val
        print("\nStep 2: Compute the product of these zeta values.")
        print(f"Product = {' * '.join(zeta_values_str)} = {zeta_product}")

    # Final calculation step
    prefactor = q / (q - 1)
    total_mass = prefactor * zeta_product

    print("\nStep 3: Combine all parts to find the total mass.")
    print(f"Total Mass = ({q} / ({q} - 1)) * {zeta_product}")
    print(f"             = {prefactor} * {zeta_product}")
    print(f"             = {total_mass}")

    # Provide exact fractional result for the example case
    if q == 2 and n == 2:
        print("\nNote: For the case q=2, n=2, the exact fractional value is 16/3.")

if __name__ == "__main__":
    main()
