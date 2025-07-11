import sympy

def solve_feynman_vertices():
    """
    Calculates the minimum number of vertices for a two-loop Feynman diagram
    in a scalar phi^4 theory.
    """

    # Define variables as symbols for symbolic manipulation
    L, I, V, E, n = sympy.symbols('L I V E n')

    # 1. The topological loop formula for a connected diagram
    loop_formula = sympy.Eq(L, I - V + 1)
    print(f"The loop formula is: {loop_formula}")

    # 2. The formula relating lines and vertices
    vertex_formula = sympy.Eq(n * V, 2 * I + E)
    print(f"The vertex formula is: {vertex_formula}")
    print("-" * 30)

    # 3. Substitute the known values for our specific problem.
    # L=2 for a two-loop diagram.
    # n=4 for a scalar phi^4 theory.
    L_val = 2
    n_val = 4
    print(f"For our problem:")
    print(f"Number of loops (L) = {L_val}")
    print(f"Lines per vertex (n) for phi^4 theory = {n_val}\n")

    # Isolate I from the loop formula
    I_solved = sympy.solve(loop_formula, I)[0]
    # This gives: I = L + V - 1

    # Substitute I into the vertex formula
    formula_with_I_substituted = vertex_formula.subs(I, I_solved)
    # This gives: n*V = 2*(L + V - 1) + E

    # Solve the resulting equation for V
    V_solved = sympy.solve(formula_with_I_substituted, V)[0]
    print(f"Solving for V, we get: V = {V_solved}")

    # Now substitute the numerical values for L and n
    V_in_terms_of_E = V_solved.subs({L: L_val, n: n_val})
    print(f"Substituting L={L_val} and n={n_val}, we have V = {V_in_terms_of_E}")
    print("-" * 30)

    # 4. To minimize V, we must minimize E. The minimum possible value for E is 0.
    E_min = 0
    print(f"To find the minimum number of vertices (V), we must use the minimum")
    print(f"possible number of external lines (E), which is E = {E_min}.")

    # 5. Calculate the minimum V
    min_V = V_in_terms_of_E.subs(E, E_min)

    # 6. Display the final calculation and the result
    print("\nFinal Calculation:")
    final_eq = f"Minimum V = (2 * (L - 1) + E) / (n - 2)"
    final_calc = f"Minimum V = (2 * ({L_val} - 1) + {E_min}) / ({n_val} - 2) = {int(min_V)}"
    print(final_eq)
    print(final_calc)
    
    print("\nThis result corresponds to a 'figure-8' vacuum diagram with one vertex,")
    print("two internal lines, and two loops.")

    # Return the final integer answer for the "answer" block
    return int(min_V)


if __name__ == '__main__':
    result = solve_feynman_vertices()
    print(f"\nThe minimum number of vertices is {result}.")
    print(f"\n<<<{result}>>>")