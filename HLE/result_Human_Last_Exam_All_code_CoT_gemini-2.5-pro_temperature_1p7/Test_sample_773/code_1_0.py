import math

def dedekind_zeta_FqT(i, q):
    """
    Calculates the value of the Dedekind zeta function Z(i) for the
    global function field F_q(T).
    Z(s) = 1 / ((1 - q^{-s}) * (1 - q^{1-s}))
    """
    if i <= 0:
        return math.nan
    if i == 1:
        return math.inf
    
    term1 = 1 - q**(-i)
    term2 = 1 - q**(1-i)
    
    if term1 == 0 or term2 == 0:
        return math.inf
        
    return 1 / (term1 * term2)

def calculate_total_mass(n, q, qv):
    """
    Calculates the total mass based on the derived formula:
    Mass = (qv / (qv - 1)) * product_{i=2 to n} Z(i)
    """
    # The product from i=2 to n is 1 if n < 2.
    z_product = 1.0
    if n >= 2:
        for i in range(2, n + 1):
            z_product *= dedekind_zeta_FqT(i, q)
            
    prefactor = qv / (qv - 1)
    total_mass = prefactor * z_product
    return total_mass

def main():
    """
    Main function to execute the calculation and print the results.
    """
    # As the problem does not specify the parameters, we use standard
    # assumptions and example values for demonstration.
    # The local field is assumed to be the completion of k=F_q(T) at infinity, so q_v = q.
    # We choose n=2 and q=3 as a non-trivial example.
    n = 2
    q = 3
    qv = q

    print(f"This script calculates the total mass for n={n}, q={q}, and q_v={qv}.")
    print("-" * 50)
    
    # Explain the formula structure as per the instructions
    print("The final calculation follows the formula:")
    
    z_prod_str = " * ".join([f"Z({i})" for i in range(2, n + 1)]) if n >= 2 else "1"
    qv_minus_1 = qv - 1
    
    print(f"Mass = ({qv} / ({qv} - 1)) * {z_prod_str}")
    print(f"     = ({qv} / {qv_minus_1}) * {z_prod_str}\n")
    
    # Show the computation for each term in the formula
    print("Let's compute the values:")
    
    z_values_str = []
    if n >= 2:
        for i in range(2, n + 1):
            val = dedekind_zeta_FqT(i, q)
            print(f"Z({i}) for q={q} is 1 / ((1 - {q}^-{i}) * (1 - {q}^(1-{i}))) = {val}")
            z_values_str.append(str(val))

    # Calculate final result
    final_mass = calculate_total_mass(n, q, qv)
    
    print("\nPutting it all together:")
    if n >= 2:
        z_prod_val = math.prod(dedekind_zeta_FqT(i, q) for i in range(2, n+1))
        print(f"Mass = ({qv} / {qv_minus_1}) * ({' * '.join(z_values_str)})")
        print(f"     = {qv / qv_minus_1} * {z_prod_val}")
        print(f"     = {final_mass}\n")
    else: # Case n=1
        print(f"Mass = {qv / qv_minus_1} * 1 = {final_mass}\n")

    print("Final answer:")
    print(final_mass)
    

if __name__ == "__main__":
    main()
