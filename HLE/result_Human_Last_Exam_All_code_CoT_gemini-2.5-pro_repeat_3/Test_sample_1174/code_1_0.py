import sys

def solve_physics_scaling():
    """
    This script calculates the value of a weighted sum of scaling exponents 
    from a physics problem concerning thermal magnetic noise.

    The exponents are determined based on physical principles and the given scaling laws:
    - n1, n2, n3 from the low-frequency limit scaling: S_B ~ sigma^n1 * T^n2 * z^n3
    - n4, n5, n6 from the frequency-dependent scaling of S_B.
    """

    # Exponents for the low-frequency limit: S_B ~ sigma^n1 * T^n2 * z^n3
    # n1: Proportionality to conductivity sigma. From fluctuation-dissipation, S_B ~ sigma.
    n1 = 1
    # n2: Proportionality to temperature T. From fluctuation-dissipation, S_B ~ T.
    n2 = 1
    # n3: Proportionality to distance z. Based on dimensional analysis and the given
    # t-independence, S_B ~ 1/z.
    n3 = -1

    # Exponents for the frequency-dependent scaling S_B(omega)
    # n4: For omega << 1/(sigma*z*t). Low-frequency white noise limit. S_B ~ const.
    n4 = 0
    # n5: For 1/(sigma*z*t) << omega << 1/(sigma*t^2). Inductive eddy current regime. S_B ~ omega^-2.
    n5 = -2
    # n6: For omega >> 1/(sigma*t^2). Skin effect regime (half-space). S_B ~ omega^-1/2.
    n6 = -0.5

    exponents = [n1, n2, n3, n4, n5, n6]
    
    # Calculate the sum: sum_{k=1 to 6} k * n_k
    total_sum = 0
    equation_parts = []
    for k in range(1, 7):
        n_k = exponents[k-1]
        term = k * n_k
        total_sum += term
        
        # Format the numbers for the equation string
        n_k_str = f"({n_k})" if n_k < 0 else f"{n_k}"
        equation_parts.append(f"{k}*{n_k_str}")

    equation_str = " + ".join(equation_parts)

    print("The exponents are:")
    print(f"n1 = {n1}")
    print(f"n2 = {n2}")
    print(f"n3 = {n3}")
    print(f"n4 = {n4}")
    print(f"n5 = {n5}")
    print(f"n6 = {n6}")
    print("\nThe calculation is for the sum of k*n_k for k from 1 to 6.")
    print(f"Sum = {equation_str}")
    print(f"Sum = {1*n1} + {2*n2} + {3*n3} + {4*n4} + {5*n5} + {6*n6}")
    print(f"Sum = {total_sum}")

solve_physics_scaling()
print("\n<<< -13 >>>")