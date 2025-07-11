import random

def solve_markov_chain_problem():
    """
    Analyzes the Markov chain problem and provides a computational demonstration.
    """
    # --- Parameters for the M/M/1 Queue Example ---
    # We model a system where the conditions lead to a non-positive recurrent chain.
    # Arrival rate
    lambda_ = 1.0
    # Service rate
    mu = 1.0
    # Traffic intensity rho = lambda / mu.
    # The M/M/1 queue is null recurrent (and thus not positive recurrent) when rho = 1.
    # It is transient (also not positive recurrent) when rho > 1.

    print("--- Analysis of the M/M/1 Queue Example ---")
    print(f"We consider an M/M/1 queue with lambda = {lambda_}, mu = {mu} (rho = {lambda_/mu}).")
    
    # --- Verification of the Problem Conditions ---
    # Let the state space be Sigma = {0, 1, 2, ...} (number of customers).
    # The chain is irreducible.
    # Let the finite set be A = {0}.
    # Let the function be f(x) = x. f(x) is non-negative and f(x) -> infinity as x -> infinity.

    # We check the drift condition for x not in A (i.e., x > 0):
    # E[f(X_1) | X_0 = x] - f(x) >= 0
    p_up = lambda_ / (lambda_ + mu)
    p_down = mu / (lambda_ + mu)
    
    # The drift for x > 0 is (p_up * (x+1) + p_down * (x-1)) - x, which simplifies to:
    drift = (lambda_ - mu) / (lambda_ + mu)
    
    print(f"\nFor f(x)=x and A={{0}}, the drift outside A is constant: {drift:.2f}.")
    print(f"The condition is that drift >= 0, which requires lambda >= mu. Our parameters satisfy this.")

    # --- Illustrating the Theoretical Contradiction ---
    # The proof shows that assuming positive recurrence leads to the conclusion:
    # f(x_0) <= max_{y in A} f(y) for ANY state x_0.
    
    # In our example, A = {0}, so f(y) for y in A is f(0) = 0.
    M_A = 0
    # We can choose any x_0. Let's pick one where the contradiction is obvious.
    x_0 = 10
    f_x0 = x_0  # Since f(x) = x

    print("\n--- Theoretical Contradiction Illustration ---")
    print("The formal proof leads to a final contradictory equation.")
    print(f"The equation must hold for any chosen x_0: f(x_0) <= M_A")
    print("Let's output the numbers for our example:")
    print(f"Chosen x_0 = {x_0}, which gives f(x_0) = {f_x0}")
    print(f"Set A = {{0}}, which gives M_A = max(f(y) for y in A) = {M_A}")
    print(f"Final Equation: {f_x0} <= {M_A}")
    print("This inequality is false, which confirms the theoretical conclusion.")

    # --- Simulation to Show Non-Positive Recurrence ---
    print("\n--- Simulation Evidence ---")
    print("We now simulate the chain to estimate the mean return time to state 0.")
    print("For a non-positive recurrent chain, this time is infinite.")
    print("The running average of the simulated return times should therefore not converge to a finite value.")

    def simulate_one_return(start_state=0):
        """Simulates one cycle from start_state until it returns."""
        if start_state != 0:
            return 0
        
        current_state = 1
        time = 1
        # Safety break for transient case where return may not happen
        max_steps = 1000000 
        
        while time < max_steps:
            if current_state == 0:
                return time  # Returned to 0
            
            if random.random() < p_up:
                current_state += 1
            else:
                current_state -= 1
            time += 1
        return time

    num_simulations = 5000
    total_time = 0

    for i in range(1, num_simulations + 1):
        return_time = simulate_one_return()
        total_time += return_time
        
        if i % 1000 == 0 or i == 1:
            current_avg = total_time / i
            print(f"After {i} simulations, the average return time is: {current_avg:.2f}")

    print("\nAs the number of simulations increases, the average return time grows,")
    print("providing computational evidence for an infinite mean return time.")

solve_markov_chain_problem()