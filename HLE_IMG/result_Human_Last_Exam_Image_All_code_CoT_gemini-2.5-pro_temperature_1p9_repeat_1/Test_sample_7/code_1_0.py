import sympy

def solve_chemistry_problem():
    """
    This function determines and prints the formula for the second voltage plateau
    of a graphite anode during lithium intercalation.

    The problem describes two voltage plateaus:
    1.  ~0.09V (between 50% and 100% stoichiometry): This corresponds to the
        transition from Stage 2 (LiC12) to Stage 1 (LiC6). The voltage is
        determined by the chemical potential for this process, mu_1.
        So, V_1 = -mu_1 / e.

    2.  ~0.13V (between 20% and 50% stoichiometry): This corresponds to the
        transition from a dilute phase to Stage 2 (LiC12). The voltage is
        determined by the chemical potential for this process, mu_2.
        So, V_2 = -mu_2 / e.

    The question asks for the formula for the second plateau.
    """

    # Define the symbols for the formula
    V = sympy.Symbol("V")
    mu = sympy.Symbol("μ")
    e = sympy.Symbol("e")
    stage_number = 2

    # The formula for the voltage of the plateau corresponding to stage k is V = -μ_k / e.
    # For the second plateau, k=2.
    formula = sympy.Eq(V, -mu/e)
    
    # Per the instruction "output each number in the final equation!", we print the
    # formula with the specific stage number.
    print(f"V = -{mu.name}_{stage_number} / {e.name}")

solve_chemistry_problem()