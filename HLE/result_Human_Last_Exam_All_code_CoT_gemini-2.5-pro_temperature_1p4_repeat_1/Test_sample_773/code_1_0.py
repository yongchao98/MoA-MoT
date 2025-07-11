import math

def calculate_mass(n, q):
    """
    Calculates the total mass based on the derived formula.
    
    The formula for the total mass is:
    Mass = (q / (q - 1)) * Product_{i=2 to n} zeta_R(i)
    where zeta_R(i) = 1 / (1 - q**(1-i)).
    """
    if n < 1 or q < 2:
        raise ValueError("n must be >= 1 and q must be >= 2.")

    if n == 1:
        # The product term is empty and equals 1.
        mass = q / (q - 1)
        print(f"For n=1, q={q}:")
        print(f"The formula is q / (q - 1)")
        print(f"Mass = {q} / ({q} - 1) = {mass}")
        return mass

    # Calculate zeta_R(i) for i from 2 to n
    zeta_values = []
    for i in range(2, n + 1):
        zeta_val = 1 / (1 - q**(1 - i))
        zeta_values.append(zeta_val)
    
    # Calculate the product of zeta values
    product_zeta = 1
    for val in zeta_values:
        product_zeta *= val

    # Prefactor q / (q - 1)
    prefactor = q / (q - 1)

    # Final mass
    mass = prefactor * product_zeta
    
    # Print the equation components
    print(f"For n={n}, q={q}:")
    print("The formula is (q / (q - 1)) * Z(2) * Z(3) * ... * Z(n)")
    
    equation_str = f"Mass = ({q} / ({q} - 1))"
    for i in range(2, n + 1):
        zeta_val = 1 / (1 - q**(1 - i))
        equation_str += f" * {zeta_val:.4f}"
    
    equation_str += f" = {mass:.4f}"
    
    print("\nCalculation:")
    print(f"Prefactor (q / (q-1)) = {prefactor:.4f}")
    for i in range(2, n + 1):
        zeta_val = 1 / (1 - q**(1 - i))
        print(f"Z({i}) = {zeta_val:.4f}")
    print(f"Product of Z(i) values = {product_zeta:.4f}")
    
    print("\nFinal Equation:")
    final_eq_parts = [f"{prefactor:.4f}"]
    for z in zeta_values:
        final_eq_parts.append(f"{z:.4f}")
    
    print(" * ".join(final_eq_parts), "=", f"{mass:.4f}")
    
    return mass

# Example usage with n=3, q=2 as the problem parameters are not specified.
n_example = 3
q_example = 2
final_mass = calculate_mass(n_example, q_example)

# The final answer in the required format would be the numerical result 
# for a specific choice of n and q. For n=3, q=2, the result is 16/3.
final_answer = 16/3