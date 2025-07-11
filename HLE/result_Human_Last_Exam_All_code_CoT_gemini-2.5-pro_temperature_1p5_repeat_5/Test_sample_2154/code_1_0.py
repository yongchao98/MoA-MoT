import math

def solve_picard_fuchs_order():
    """
    Calculates the minimal order u_r(n) of the Picard-Fuchs differential
    equation for the Hamiltonian V(q) = 1/2 * (q^2 - q^n) for n = 3 to 12.
    """
    final_results = []
    print("This script calculates the minimal order u_r(n) of the Picard-Fuchs differential equation for the given Hamiltonian.")
    print("The formula for the order depends on the symmetry of the potential V(q), which is determined by the parity of n.")
    print("-" * 60)

    for n in range(3, 13):
        if n % 2 != 0:  # n is odd
            print(f"Calculating for n = {n} (odd case):")
            print("  The potential V(q) is not symmetric, so we use the general formula.")
            
            g_numerator = n - 1
            g_denominator = 2
            g = g_numerator // g_denominator
            order = 2 * g
            
            print(f"  The genus of the associated hyperelliptic curve is g = floor((n-1)/2).")
            print(f"  g = floor(({n}-1)/2) = floor({g_numerator}/{g_denominator}) = {g}.")
            print(f"  The order of the Picard-Fuchs equation is u_r(n) = 2 * g.")
            print(f"  u_r({n}) = 2 * {g} = {order}")
            final_results.append(order)
        else:  # n is even
            print(f"Calculating for n = {n} (even case):")
            print("  The potential V(q) is symmetric, which leads to a reduction in the order of the equation.")

            numerator = n - 2
            denominator = 4
            intermediate_val = numerator / denominator
            order = math.ceil(intermediate_val)
            
            print(f"  The order is given by the reduced formula u_r(n) = ceil((n-2)/4).")
            print(f"  u_r({n}) = ceil(({n}-2)/4) = ceil({numerator}/{denominator}) = ceil({intermediate_val}) = {order}")
            final_results.append(order)
            
        print("-" * 60)

    print("The complete set of values is:")
    print(f"{{u_r(3), u_r(4), ..., u_r(12)}} = {final_results}")


solve_picard_fuchs_order()