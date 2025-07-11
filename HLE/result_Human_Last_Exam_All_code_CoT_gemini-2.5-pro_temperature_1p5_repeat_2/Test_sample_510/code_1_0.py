import sympy
from sympy import symbols, Integral, Rational, log, sqrt, pretty_print

def solve_queueing_problem():
    """
    Calculates the optimal mean response time 'x' for the given M/G/1 SRPT queue.
    """
    # Define symbols and constants
    s, y = symbols('s y', real=True, positive=True)
    lambda_val = Rational(3, 2)
    job_size_upper_bound = 1

    # PDF of job size is f(s) = 1 for s in [0, 1]
    f_s = 1

    # 1. Calculate the mean service time E[S]
    E_S = Integral(s * f_s, (s, 0, job_size_upper_bound)).doit()

    # 2. Calculate rho(s), the load from jobs smaller than s
    # rho(s) = lambda * integral from 0 to s of y*f(y) dy
    rho_s = lambda_val * Integral(y * f_s, (y, 0, s)).doit()

    # 3. Calculate the mean waiting time E[W]
    # The integrand is s^2*f(s) / (1 - rho(s))^2
    integrand_W = (s**2 * f_s) / (1 - rho_s)**2
    # E[W] = (lambda/2) * integral from 0 to 1 of the integrand
    E_W = (lambda_val / 2) * Integral(integrand_W, (s, 0, job_size_upper_bound)).doit()

    # 4. Calculate the optimal mean response time x = E[T] = E[W] + E[S]
    x = E_W + E_S

    # Print the results step-by-step
    print("Step 1: Calculate the mean job size E[S]")
    print(f"E[S] = integral from 0 to 1 of (s * 1) ds")
    print(f"E[S] = {E_S}")
    print("-" * 30)

    print("Step 2: Calculate the mean waiting time E[W]")
    print(f"Using lambda = {lambda_val}")
    print(f"rho(s) = {lambda_val} * integral from 0 to s of (y * 1) dy = {rho_s}")
    print(f"E[W] = ({lambda_val}/2) * integral from 0 to 1 of [s^2 / (1 - {rho_s})^2] ds")
    
    # Use pretty_print for better formatting of complex expressions
    print("\nValue of E[W]:")
    pretty_print(E_W)
    print("-" * 30)

    print("Step 3: Calculate the optimal mean response time x = E[T] = E[W] + E[S]")
    print(f"x = ({E_W}) + ({E_S})")

    print("\nFinal expression for x:")
    pretty_print(x)
    print("\nFormatted equation:")
    # Using SymPy's log as ln for standard notation
    print(f"x = {x.args[0]} + {x.args[1].args[0]}*{x.args[1].args[1]}*ln({x.args[1].args[2].args[0]} + {x.args[1].args[2].args[1]})")

solve_queueing_problem()