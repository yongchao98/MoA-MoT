import sympy

def calculate_final_mass():
    """
    Calculates the total mass requested in the problem statement.
    
    This function prompts the user for the necessary parameters and then
    computes the result based on established formulas for volumes of 
    S-arithmetic quotients in function fields.
    """
    try:
        n_str = input("Enter the dimension n (integer >= 1): ")
        n = int(n_str)
        if n < 1:
            raise ValueError("n must be a positive integer.")

        q_str = input("Enter q, the size of the constant field of R (prime power > 1): ")
        q = int(q_str)
        if not sympy.isprime(q) and not sympy.is_prime_power(q):
            print("Warning: q is typically a prime power.")
        if q <= 1:
            raise ValueError("q must be greater than 1.")
            
        q_v_str = input(f"Enter q_v, the order of the residual field of K_hat (a power of {q}): ")
        q_v = int(q_v_str)
        if q_v <= 1:
            raise ValueError("q_v must be greater than 1.")

        P_K_coeffs_str = input("Enter the coefficients of the zeta polynomial P_K(u) as a comma-separated list (e.g., '1' for g=0, or '1,2,3' for 1+2u+3u^2): ")
        P_K_coeffs = [int(c.strip()) for c in P_K_coeffs_str.split(',')]
        
    except (ValueError, TypeError) as e:
        print(f"Invalid input. Please enter valid parameters. Error: {e}")
        return

    u = sympy.Symbol('u')
    P_K = sum(c * u**i for i, c in enumerate(P_K_coeffs))

    print("\n--- Calculation Steps ---")
    print(f"Parameters: n = {n}, q = {q}, q_v = {q_v}, P_K(u) = {P_K}")

    if n == 1:
        # For n=1, the product of zeta functions is empty, conventionally 1.
        mu_X = sympy.Integer(1)
        print("\nFor n=1, the volume of the space of lattices mu(X) is 1.")
        zeta_values_for_print = []
    else:
        def calculate_zeta_K(s, q_val, poly):
            """Calculates the Dedekind zeta function zeta_K(s)."""
            q_inv_s = sympy.Pow(q_val, -s)
            numerator = poly.subs(u, q_inv_s)
            denominator = (1 - q_inv_s) * (1 - q_val * q_inv_s)
            return numerator / denominator

        print("\nCalculating values of the Dedekind zeta function Z(k) = zeta_K(k):")
        zeta_values = []
        for k in range(2, n + 1):
            val = calculate_zeta_K(k, q, P_K)
            zeta_values.append(val)
            print(f"Z({k}) = {val.evalf()} (Exact: {val})")

        mu_X = sympy.prod(zeta_values)
        zeta_values_for_print = zeta_values
        print(f"\nThe volume mu(X) = Z(2) * ... * Z({n}) is {mu_X.evalf()} (Exact: {mu_X})")
    
    factor = (sympy.Integer(q_v) * (sympy.Integer(q) - 1)) / (sympy.Integer(q_v) - 1)
    print(f"\nThe multiplier factor is {q_v}({q}-1)/({q_v}-1) = {factor.evalf()} (Exact: {factor})")
    
    total_mass = factor * mu_X
    
    print("\n--- Final Answer ---")
    print(f"The total mass is the multiplier multiplied by the volume mu(X).")

    # Format the final equation as requested
    factor_str = f"({q_v}*({q}-1)/({q_v}-1))"
    mass_val_str = f"{mu_X.evalf()}"
    if isinstance(mu_X, sympy.Expr) and not mu_X.is_Number:
       mass_val_str += f" (Exact: {mu_X})"
       
    if n > 1:
        zeta_vals_str = " * ".join([f"Z({k})" for k in range(2, n + 1)])
        print(f"Total Mass = {factor_str} * [{zeta_vals_str}]")
        
        zeta_vals_num_str = " * ".join([f"{z.evalf()}" for z in zeta_values_for_print])
        print(f"             = {factor.evalf()} * [{zeta_vals_num_str}]")
        print(f"             = {factor.evalf()} * {mu_X.evalf()}")
    else:
        print(f"Total Mass = {factor_str} * 1")
    
    print(f"             = {total_mass.evalf()}")

    # Finally, output the final numeric equation as required by the prompt instructions
    print("\nFinal numeric equation:")
    if n > 1:
        zeta_vals_num_str = " * ".join([f"{z.evalf()}" for z in zeta_values_for_print])
        print(f"{factor.evalf()} * {zeta_vals_num_str} = {total_mass.evalf()}")
    else:
        print(f"{factor.evalf()} * 1 = {total_mass.evalf()}")

if __name__ == '__main__':
    calculate_final_mass()