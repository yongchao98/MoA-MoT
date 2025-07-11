import math

def calculate_prefactor_cn():
    """
    This script determines the prefactor c_n for the fully f-connected diagram's
    contribution to the n-th virial coefficient B_n.

    The n-th virial coefficient B_n is given by the Mayer cluster theory as:
    
        B_n = - (n-1)/n * beta_{n-1}
    
    where beta_{n-1} is the (n-1)-th irreducible cluster integral. It is defined as:
    
        beta_{n-1} = (1/(n-1)!) * Integral(d(r_2)...d(r_n)) * [Sum over all biconnected graphs on n labeled vertices of Product(f_ij)]
    
    Combining these gives:
    
        B_n = - (n-1)/n! * Integral(d(r)^(n-1)) * [Sum over biconnected graphs G_n of Product(f_ij)]
    
    The fully f-connected diagram is the complete graph K_n, which is one of the biconnected graphs in the sum. Its contribution to B_n is:
    
        Contribution(K_n) = - (n-1)/n! * Integral(d(r)^(n-1)) * [Product_{i<j} f_ij]
    
    The user defines this contribution as c_n * Lambda_n, where:
    
        Lambda_n = Integral(d(r)^(n-1)) * [Product_{i<j} f_ij]
        
    By comparing the two expressions, we can identify the prefactor c_n.
    """
    
    print("Derivation of the prefactor c_n:")
    print("The contribution of the fully f-connected graph (K_n) to the n-th virial coefficient (B_n) is given by the term:")
    print("    Contribution = -((n-1)/n!) * Integral( Product_{i<j} f_ij )")
    print("Comparing this with the definition B_n = c_n * Lambda_n + ..., where Lambda_n is the integral part, we get:")
    print("    c_n = -(n - 1) / n!")
    print("\nCalculating c_n for the first few values of n:")
    print("--------------------------------------------------")

    for n in range(2, 6):
        if n < 2:
            print(f"For n = {n}: B_n is not defined.")
            continue

        numerator = -(n - 1)
        # In Python, math.factorial calculates n!
        denominator = math.factorial(n)
        
        # To print a simplified fraction, we find the greatest common divisor.
        common_divisor = math.gcd(abs(numerator), denominator)
        simplified_num = numerator // common_divisor
        simplified_den = denominator // common_divisor

        value = numerator / denominator

        print(f"For n = {n}:")
        # Outputting each number in the final equation, as requested.
        print(f"  c_{n} = -({n} - 1) / {n}!")
        print(f"      = {numerator} / {denominator}")
        # Only print the simplified fraction if it's different.
        if denominator != simplified_den:
             print(f"      = {simplified_num} / {simplified_den}")
        print(f"      = {value}")
        print("-" * 20)

# Execute the function to perform the calculation and print the results.
calculate_prefactor_cn()