import sympy

def solve_steady_state():
    """
    This function derives the steady-state probability pi_0 for the given
    birth-death process and prints the derivation steps.
    """

    # Define rho as a symbolic variable for clear output
    rho = sympy.Symbol('rho')
    pi_0 = sympy.Symbol('pi_0')
    i = sympy.Symbol('i')
    k = sympy.Symbol('k')

    print("Step 1: State the detailed balance equations for a birth-death process.")
    print("pi_i * lambda_i = pi_{i+1} * mu_{i+1}\n")

    print("Step 2: Substitute the given rates into the equations.")
    print("Arrival rate: lambda_i = lambda / (i+1)")
    print("Service rate: mu_i = mu")
    print(f"So, for state i+1, mu_{'i+1'} = mu")
    print("The equation becomes: pi_i * (lambda / (i+1)) = pi_{i+1} * mu\n")

    print(f"Step 3: Express pi_{'i+1'} in terms of pi_i and rho = lambda/mu.")
    print("pi_{i+1} = pi_i * (lambda / (mu * (i+1)))")
    print(f"pi_{'i+1'} = pi_i * ({rho} / (i+1))\n")

    print("Step 4: Find a general expression for pi_k in terms of pi_0.")
    print(f"For i=0: pi_1 = pi_0 * ({rho} / 1)")
    print(f"For i=1: pi_2 = pi_1 * ({rho} / 2) = (pi_0 * {rho}) * ({rho} / 2) = pi_0 * {rho}**2 / 2!")
    print(f"For i=2: pi_3 = pi_2 * ({rho} / 3) = (pi_0 * {rho}**2 / 2) * ({rho} / 3) = pi_0 * {rho}**3 / 3!")
    print("By induction, the general formula is:")
    print(f"pi_k = pi_0 * ({rho}**k / k!)\n")

    print("Step 5: Use the normalization condition Sum(pi_k) = 1 for k=0 to infinity.")
    print("Sum(pi_0 * (rho**k / k!)) = 1")
    print("pi_0 * Sum(rho**k / k!) = 1\n")

    print("Step 6: Recognize that Sum(rho**k / k!) is the Taylor series for e**rho.")
    print("The sum from k=0 to infinity of (x**k / k!) is e**x.")
    print(f"So, Sum({rho}**k / k!) = e**{rho}\n")

    print("Step 7: Substitute this back and solve for pi_0.")
    print(f"pi_0 * e**{rho} = 1")
    print(f"pi_0 = 1 / e**{rho}")
    print("\n--- Final Answer ---")
    print("The steady-state probability pi_0 is:")
    print(f"pi_0 = e**(-{rho})")


solve_steady_state()