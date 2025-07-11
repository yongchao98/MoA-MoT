import math

def calculate_total_mass(n, q):
    """
    Calculates the total mass for the given n and q based on the problem statement.

    Args:
        n (int): The dimension of the general linear group, must be a positive integer.
        q (int): The characteristic of the local field, must be an integer > 1.
    """
    # We assume the most common setting where the residue field size q_v equals the characteristic q.
    q_v = q
    
    print(f"Calculating the total mass for n = {n}, q = {q}, assuming q_v = q.")

    # --- Step 1: Define the formula components ---
    # The Dedekind zeta function for the rational function field F_q(T) is Z(s).
    def Z(s, q_val):
        return 1.0 / ((1.0 - q_val**-s) * (1.0 - q_val**(1.0-s)))

    # The volume of the space of lattices mu(X) is given by:
    # mu(X) = (1 / (q - 1)) * Product_{k=2 to n} Z(k)
    # The total mass M is (q_v * (q - 1)) / (q_v - 1) * mu(X)
    # With q_v = q, this simplifies to M = (q / (q - 1)) * Product_{k=2 to n} Z(k)
    
    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.")
        return
    if not isinstance(q, int) or q < 2:
        print("Error: q must be an integer greater than 1.")
        return

    # --- Step 2: Calculate the product of zeta values ---
    product_Z_k = 1.0
    if n >= 2:
        for k in range(2, n + 1):
            product_Z_k *= Z(k, q)

    # --- Step 3: Calculate the final total mass ---
    # The coefficient is q_v / (q_v - 1)
    coefficient = q_v / (q_v - 1.0)
    total_mass = coefficient * product_Z_k

    # --- Step 4: Output the equation and its components as requested ---
    print("\nThe equation for the total mass is:")
    
    # Build and print the equation string with symbols
    equation_str = f"Mass = ({q_v} / ({q_v} - 1))"
    if n >= 2:
        equation_str += " * " + " * ".join([f"Z({k})" for k in range(2, n + 1)])
    print(equation_str)

    # Print the values of each number in the equation
    print("\nComponent values:")
    print(f"  Coefficient = {q_v} / ({q_v} - 1) = {coefficient}")
    if n >= 2:
        print("  Zeta function values:")
        for k in range(2, n + 1):
            print(f"    Z({k}) = {Z(k, q)}")
        print(f"  Product of Zeta values = {product_Z_k}")
    
    print(f"\nFinal Resulting Total Mass = {total_mass}")

# You can modify these values to compute the mass for different n and q.
n_val = 3
q_val = 2

calculate_total_mass(n_val, q_val)