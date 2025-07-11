import sympy
from sympy import S, Symbol, limit, oo, Eq, solve

def find_infimum_c():
    """
    This function calculates the infimum of c for which the given Markov chain is transient.
    It uses sympy for symbolic calculations to apply the transience criterion for a 1D random walk.
    """
    # Define symbolic variables for the state k and the parameter c
    k = Symbol('k', positive=True, integer=True)
    c = Symbol('c')

    # Define the possible jumps and their corresponding probabilities for large k
    jumps = [-2, -1, 1, 2]
    probs = [S(1)/4, S(1)/4 - c/k, S(1)/4 + c/k, S(1)/4]

    # Step 1: Calculate the expected drift (mean jump), mu_k.
    # For a drift mu_k ~ A/k, we find the coefficient A.
    mu_k = sympy.simplify(sum(j * p for j, p in zip(jumps, probs)))
    A = sympy.simplify(mu_k * k)

    # Step 2: Calculate the second moment of the jump and its limit B as k -> infinity.
    second_moment_k = sympy.simplify(sum(j**2 * p for j, p in zip(jumps, probs)))
    B = limit(second_moment_k, k, oo)

    print("For a Markov chain with drift mu_k ~ A/k and limiting jump variance B,")
    print("the chain is transient if A > B/2.")
    print("The critical condition separating transience from recurrence is A = B/2.")
    
    print("\nBased on the chain's transition probabilities, we calculate:")
    print(f"Coefficient of drift, A = {A}")
    print(f"Limiting variance, B = {B}")

    # Step 3: Set up and solve the critical equation A = B / 2.
    critical_value = B / 2
    
    print("\nThe critical equation is A = B/2, which numerically is:")
    # The prompt requires outputting each number in the final equation.
    # The symbolic expression for A is printed, and the numeric values for B and B/2.
    print(f"{A} = {B} / 2")
    print("This simplifies to the final equation:")
    print(f"{A} = {critical_value}")
    
    # Solve the equation for c to find the infimum.
    critical_equation = Eq(A, critical_value)
    solution = solve(critical_equation, c)
    infimum_c = solution[0]

    print(f"\nSolving for c gives the infimum value:")
    print(f"c = {float(infimum_c)}")

find_infimum_c()

print("<<<0.625>>>")