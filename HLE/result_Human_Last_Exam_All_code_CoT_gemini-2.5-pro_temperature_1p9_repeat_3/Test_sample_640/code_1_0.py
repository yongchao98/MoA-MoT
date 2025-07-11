def solve_scattering_statements():
    """
    Analyzes the logical statements about quantum scattering matrices and
    prints the numbers of the correct statements in a final equation.

    The logic is based on established theorems in multichannel scattering theory:
    1. The coupling nature (trivial/nontrivial) of the potential V(r) and the
       Jost matrix F(E) are equivalent due to the uniqueness of the inverse
       scattering problem (F -> V).
       So: V is trivially coupled <=> F is trivially coupled.

    2. A trivially coupled Jost matrix F(E) algebraically yields a trivially
       coupled S-matrix S(E).
       So: F is trivially coupled => S is trivially coupled.

    3. Nontrivially coupled "transparent" potentials exist. These have a
       nontrivially coupled V(r) but a diagonal (and thus trivially coupled)
       S-matrix. This serves as a key counterexample for some implications.
    """

    # Analysis of each statement based on the logic above.
    # Stmt 1: S_NTC => V_NTC. (Correct, as it's the contrapositive of V_TC => S_TC)
    # Stmt 2: S_diag => V_diag. (Incorrect, transparent potentials are a counterexample)
    # Stmt 3: V_NTC => F_NTC. (Correct, from V_TC <=> F_TC)
    # Stmt 4: F_NTC => S_NTC. (Incorrect, transparent potentials have F_NTC but S_TC)
    # Stmt 5: Exists V_NTC with F_diag. (Incorrect, F_diag => F_TC => V_TC, a contradiction)

    correct_statements = [1, 3]

    print("The following statements are correct:")
    for statement_number in correct_statements:
        print(f"Statement {statement_number}")

    # The prompt creatively asks to "output each number in the final equation".
    # We will construct a simple sum to satisfy this format.
    equation_string = " + ".join(map(str, correct_statements))
    total_sum = sum(correct_statements)
    
    print("\nThe final 'equation' containing the numbers of the correct statements is:")
    print(f"{equation_string} = {total_sum}")


solve_scattering_statements()