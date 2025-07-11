def solve_neron_severi_rank():
    """
    Calculates the smallest and largest possible rank of the Neron-Severi group
    of the 15th symmetric power of a genus 3 Riemann surface.
    """
    g = 3

    # --- Smallest Rank Calculation ---
    # For a generic curve, End(J(C)) = Z, which has rank 1.
    # The formula for the rank of the Neron-Severi group is rho = 1 + rank(End(J(C))).
    rank_end_min = 1
    rho_min = 1 + rank_end_min

    print("Calculation for the smallest possible rank:")
    print(f"The genus of the curve C is g = {g}.")
    print("For a generic curve, the rank of the endomorphism ring of its Jacobian, rank(End(J(C))), is 1.")
    print(f"Smallest rank = 1 + rank(End(J(C))) = 1 + {rank_end_min} = {rho_min}")
    print("-" * 20)

    # --- Largest Rank Calculation ---
    # For a special curve with maximal complex multiplication, J(C) is isogenous to E^g,
    # where E is a CM elliptic curve.
    # rank(End(E)) for a CM elliptic curve is 2.
    # rank(End(J(C))) is g^2 * rank(End(E)).
    rank_end_E = 2
    rank_end_max = g**2 * rank_end_E
    rho_max = 1 + rank_end_max

    print("Calculation for the largest possible rank:")
    print(f"The genus of the curve C is g = {g}.")
    print("For a curve with maximal CM, J(C) is isogenous to E^g, where E is a CM elliptic curve.")
    print(f"The rank of End(E) for a CM elliptic curve is {rank_end_E}.")
    print(f"The maximum rank of End(J(C)) is g^2 * rank(End(E)) = {g}^2 * {rank_end_E} = {rank_end_max}.")
    print(f"Largest rank = 1 + rank(End(J(C))) = 1 + {rank_end_max} = {rho_max}")
    print("-" * 20)
    
    print(f"The smallest possible rank is {rho_min}.")
    print(f"The largest possible rank is {rho_max}.")

if __name__ == '__main__':
    solve_neron_severi_rank()
