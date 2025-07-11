import sympy

def solve_production_problem():
    """
    Calculates the probability that the production process will reach a point 
    where exactly 50% of the products are good using a martingale approach.
    """
    # 1. Initial State
    W0, B0 = 2, 1
    
    # 2. Define the Martingale and its initial value
    # M_t = (W_t - B_t) / (W_t + B_t)
    M0 = sympy.S(W0 - B0) / sympy.S(W0 + B0)
    
    # 3. Analyze the case where the process never stops (T = infinity)
    # The proportion of good products X = W_t/(W_t+B_t) converges to a Beta distribution.
    # The PDF for Beta(a, b) is x**(a-1)*(1-x)**(b-1) / Beta(a,b).
    # For Beta(2, 1), the PDF f(x) = 2x for x in [0, 1].
    x = sympy.symbols('x')
    pdf = 2 * x
    
    # If the process never stops, W_t > B_t for all t, which means X >= 1/2.
    # We need to calculate the conditional expectation E[M_inf | T=inf]
    # M_inf = limit of M_t = limit of (W_t/N_t - B_t/N_t) = X - (1-X) = 2X - 1.
    
    # Calculate the probability of the condition X >= 1/2
    prob_condition = sympy.integrate(pdf, (x, sympy.S.One/2, 1))
    
    # The conditional PDF is pdf / prob_condition
    conditional_pdf = pdf / prob_condition
    
    # Calculate the conditional expectation of M_inf
    E_M_inf_conditional = sympy.integrate((2*x - 1) * conditional_pdf, (x, sympy.S.One/2, 1))

    # 4. Solve for the probability p = P(T < infinity)
    # The Optional Stopping Theorem gives: M0 = (1 - p) * E[M_inf | T=inf]
    # Here, p is the probability of stopping.
    p = sympy.symbols('p')
    equation = sympy.Eq(M0, (1 - p) * E_M_inf_conditional)
    solution = sympy.solve(equation, p)
    
    p_final = solution[0]
    
    # 5. Print the results step-by-step
    print("--- Martingale Solution ---")
    print(f"Initial state: Good products W_0 = {W0}, Defective products B_0 = {B0}")
    print(f"Martingale M_t = (W_t - B_t) / (W_t + B_t)")
    print(f"Initial value M_0 = ({W0} - {B0}) / ({W0} + {B0}) = {M0}")
    print("\nIf the process never stops (T = infinity):")
    print(f"The limiting proportion of good products X follows a Beta({W0}, {B0}) distribution.")
    print(f"The condition for not stopping is X >= 1/2.")
    print(f"The probability of this condition is P(X >= 1/2) = {prob_condition}")
    print(f"The expected value of the martingale M_t in this case is E[M_inf | T=inf] = {E_M_inf_conditional}")
    print("\nUsing the Optional Stopping Theorem, we form the equation:")
    print(f"M_0 = (1 - p) * E[M_inf | T=inf]")
    print(f"Substituting the values: {M0} = (1 - p) * {E_M_inf_conditional}")
    print(f"\nSolving for p, the probability of stopping:")
    print(f"p = {p_final}")

solve_production_problem()
<<<2/5>>>