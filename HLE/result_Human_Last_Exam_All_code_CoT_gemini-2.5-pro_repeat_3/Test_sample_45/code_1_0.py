import sympy as sp

def solve_telescope_problem():
    """
    This function symbolically derives the relationship between the focal length 'f'
    of a liquid-mirror telescope and time 't' when driven by a constant power source.
    It computes the exponent 'n' in the relation f ∝ t^n.
    """
    # --- Step 1: Define all symbols ---
    # t: time
    # P: constant power
    # I: moment of inertia (constant)
    # g: acceleration due to gravity
    # C1: integration constant
    # omega_squared: a symbol to represent the square of angular velocity
    t, P, I, g, C1 = sp.symbols('t P I g C1', positive=True, real=True)
    omega_squared_symbol = sp.Symbol('omega_squared')
    f_symbol = sp.Symbol('f')
    
    # Let omega be a function of time, but we will mostly work with its integrated form.
    omega = sp.Function('omega')(t)

    # --- Step 2: Establish the physical relationships as equations ---

    # Relationship 1: Focal length 'f' and angular velocity 'omega'.
    # For a rotating liquid mirror, f = g / (2 * omega^2).
    f_vs_omega_eq = sp.Eq(f_symbol, g / (2 * omega_squared_symbol))
    
    # Relationship 2: Angular velocity 'omega' and time 't' under constant power 'P'.
    # The governing differential equation is P = I * omega * d(omega)/dt.
    # Separating variables gives P * dt = I * omega * d(omega).
    # Integrating both sides yields: P*t = (1/2)*I*omega^2 + C1
    integrated_eq = sp.Eq(P * t, I * omega(t)**2 / 2 + C1)

    # --- Step 3: Solve the equations ---
    
    # Use the initial condition omega(0) = 0 (starts from rest) to find C1.
    initial_cond_eq = integrated_eq.subs({t: 0, omega(t): 0})
    C1_val = sp.solve(initial_cond_eq, C1)[0]

    # Substitute C1 back to get the equation for omega^2(t).
    omega_eq_of_t = integrated_eq.subs(C1, C1_val)

    # Solve for omega^2 as an expression of time 't'.
    omega_squared_expr = sp.solve(omega_eq_of_t, omega(t)**2)[0]

    # --- Step 4: Combine relationships to find f(t) ---
    
    # Substitute the expression for omega^2 into the focal length equation.
    f_expr_of_t = f_vs_omega_eq.rhs.subs(omega_squared_symbol, omega_squared_expr)
    
    # --- Step 5: Extract the exponent 'n' ---

    # The expression for f(t) is proportional to t^n. We extract the exponent of t.
    # The as_powers_dict() method returns a dictionary of {base: exponent}.
    powers_dict = f_expr_of_t.as_powers_dict()
    n = powers_dict.get(t, 0) # Use .get() to be safe.

    # --- Print the results and derivation steps ---
    print("Derivation of the exponent 'n' for f ∝ t^n:")
    print("-" * 50)
    print("1. The focal length 'f' of a rotating liquid mirror is related to its angular velocity 'ω' by:")
    print(f"   f = g / (2 * ω^2)")
    print("\n2. Under a constant power source 'P', the relationship between 'ω' and time 't' is found by solving P = I * ω * dω/dt.")
    print("   With the initial condition ω(0) = 0, we find:")
    print(f"   ω^2 = {omega_squared_expr}")
    print("\n3. Substituting the expression for ω^2 into the focal length equation gives 'f' as a function of 't':")
    final_f_eq = sp.Eq(f_symbol, f_expr_of_t)
    print(f"   {final_f_eq}")
    print("\n4. The final equation shows that f is proportional to 1/t, which can be written as t^-1.")
    print("   Therefore, the exponent n in the relation f ∝ t^n is:")
    print(f"   n = {int(n)}")

solve_telescope_problem()