import sympy
from sympy import symbols, Eq, solve, Sum, oo, factorial, exp

def solve_steady_state_probability():
    """
    This function calculates the steady-state probability pi_0 for a given
    birth-death process using symbolic mathematics.
    """
    # Define the symbolic variables
    # pi_0 is the steady-state probability of being in state 0
    # rho is the ratio of arrival rate to service rate (lambda/mu)
    # k is the state index
    pi_0 = symbols('pi_0')
    rho = symbols('rho', positive=True)
    k = symbols('k', integer=True, nonnegative=True)

    # From the detailed balance equations for this process, we derive the
    # general expression for the steady-state probability of being in state k:
    # pi_k = pi_0 * (rho^k / k!)
    # This is because pi_{i+1} = pi_i * (lambda_i / mu_{i+1})
    # pi_{i+1} = pi_i * ( (lambda/(i+1)) / mu ) = pi_i * (rho / (i+1))
    # pi_1 = pi_0 * (rho/1)
    # pi_2 = pi_1 * (rho/2) = pi_0 * (rho^2 / 2!)
    # pi_k = pi_0 * (rho^k / k!)

    pi_k_term = pi_0 * rho**k / factorial(k)

    # The normalization condition states that the sum of all probabilities must be 1.
    # Sum_{k=0 to infinity} pi_k = 1
    # Sum_{k=0 to infinity} pi_0 * (rho^k / k!) = 1
    # pi_0 * Sum_{k=0 to infinity} (rho^k / k!) = 1

    # The sum is the Taylor series expansion for e^rho.
    # So, pi_0 * e^rho = 1

    # We can use sympy to perform this summation and solve the equation.
    total_probability_sum = Sum(pi_k_term, (k, 0, oo))

    # Evaluate the sum
    evaluated_sum = total_probability_sum.doit()

    # Create the normalization equation
    normalization_equation = Eq(evaluated_sum, 1)

    # Solve for pi_0
    solution = solve(normalization_equation, pi_0)

    # The solution is a list, so we extract the first element.
    final_expression = solution[0]

    # Print the final equation for pi_0
    # The prompt asks to output each number in the final equation.
    # The final equation is pi_0 = e^(-rho). There are no explicit numbers,
    # but we can represent it as pi_0 = 1 * e^(-1 * rho).
    base = sympy.E
    coefficient = -1
    print(f"The final equation is pi_0 = {base}^({coefficient}*rho)")
    print(f"In a more standard format: pi_0 = {final_expression}")


if __name__ == '__main__':
    solve_steady_state_probability()