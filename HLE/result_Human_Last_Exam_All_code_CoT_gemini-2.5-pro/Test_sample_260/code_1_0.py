def solve_hurewicz_rank():
    """
    Calculates the rank of the kernel of the Hurewicz homomorphism for the given space Y.
    """
    # Orders of the fundamental groups of X1, X2, and X3.
    # pi_1(X1) = Z_n1, pi_1(X2) = Z_n2, pi_1(X3) = Z_n3
    n1 = 5  # From the pentagon
    n2 = 8  # From the octagon
    n3 = 2  # From the real projective plane

    print("We want to find the rank of the kernel K of the Hurewicz homomorphism h*: pi_1(Y) -> H_1(Y).")
    print("The space Y is the connected sum Y = X1 # X2 # X3.")
    print("\nStep 1: Define the fundamental groups and their orders.")
    print(f"pi_1(X1) = Z_{n1}, pi_1(X2) = Z_{n2}, pi_1(X3) = Z_{n3}.")
    print(f"pi_1(Y) is the free product Z_{n1} * Z_{n2} * Z_{n3}.")
    print(f"H_1(Y) is the direct sum Z_{n1} + Z_{n2} + Z_{n3}.\n")

    print("Step 2: Use the Euler characteristic formula for the short exact sequence 1 -> K -> pi_1(Y) -> H_1(Y) -> 1.")
    print("The formula is: chi(pi_1(Y)) = chi(K) * chi(H_1(Y))\n")

    print("Step 3: Calculate each term.")
    # For a free group K of rank r, chi(K) = 1 - r.
    print("The kernel K is a free group. If its rank is r, then chi(K) = 1 - r.")
    
    # Calculate order and chi of H_1(Y)
    order_h1 = n1 * n2 * n3
    print(f"The order of H_1(Y) is {n1} * {n2} * {n3} = {order_h1}.")
    print(f"The Euler characteristic chi(H_1(Y)) = 1 / |H_1(Y)| = 1 / {order_h1}.\n")

    # Calculate chi of pi_1(Y)
    print(f"The Euler characteristic chi(pi_1(Y)) = chi(Z_{n1}) + chi(Z_{n2}) + chi(Z_{n3}) - (3-1)")
    print(f"chi(pi_1(Y)) = 1/{n1} + 1/{n2} + 1/{n3} - 2")
    
    chi_pi1_numerator = (n2 * n3) + (n1 * n3) + (n1 * n2) - (2 * n1 * n2 * n3)
    chi_pi1_denominator = n1 * n2 * n3
    
    print(f"  = ({n2 * n3} + {n1 * n3} + {n1 * n2} - {2 * n1 * n2 * n3}) / {chi_pi1_denominator}")
    print(f"  = {chi_pi1_numerator}/{chi_pi1_denominator}\n")

    print("Step 4: Form the final equation and solve for the rank r.")
    print("The equation is: chi(pi_1(Y)) = (1 - r) * chi(H_1(Y))")
    print(f"Substituting the calculated values, we get the final equation:")
    print(f"{chi_pi1_numerator}/{chi_pi1_denominator} = (1 - r) * (1/{order_h1})")
    
    print("\nSolving for r:")
    print(f"{chi_pi1_numerator} = 1 - r")
    rank = 1 - chi_pi1_numerator
    print(f"r = 1 - ({chi_pi1_numerator})")
    print(f"r = {rank}")

if __name__ == '__main__':
    solve_hurewicz_rank()
<<<95>>>