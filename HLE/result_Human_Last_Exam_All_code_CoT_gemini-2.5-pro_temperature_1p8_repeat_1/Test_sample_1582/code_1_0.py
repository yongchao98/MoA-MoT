import random

def simulate_non_positive_recurrent_chain():
    """
    This function simulates a Markov chain that fits the described properties
    to demonstrate it is not positive recurrent.

    The Markov Chain:
    - State space Σ = {0, 1, 2, ...}.
    - Transitions:
      - From state 0, always move to 1: p(0, 1) = 1.
      - From state x > 0, move to x+1 with probability p=0.6, and to x-1 with probability q=0.4.
    - This chain is irreducible.

    The Function f and Set A:
    - Let A = {0}, which is a finite set.
    - Let f(x) = 1.5^x. This function is non-negative and f(x) -> infinity as x -> infinity.

    Checking the condition sum(p(x,y)f(y)) - f(x) >= 0 for x not in A:
    For x > 0:
    E[f(X_1) | X_0=x] - f(x) = [p * f(x+1) + q * f(x-1)] - f(x)
                            = [0.6 * (1.5)^(x+1) + 0.4 * (1.5)^(x-1)] - (1.5)^x
                            = (1.5)^x * [0.6 * 1.5 + 0.4 / 1.5 - 1]
                            = (1.5)^x * [0.9 + 0.266... - 1]
                            = (1.5)^x * [0.166...] > 0.
    The condition holds.
    
    The simulation will show that many paths starting from 0 never return,
    implying the mean return time is infinite. This is a hallmark of a chain
    that is not positive recurrent.
    """
    
    # Parameters for simulation
    num_simulations = 10000
    max_steps = 2000 # A boundary to detect escape to infinity
    
    returned_count = 0
    total_return_time = 0

    print("Starting simulation for a non-positive recurrent chain...")

    for i in range(num_simulations):
        # Start each path from state 0. By definition, it must first go to 1.
        state = 1
        time = 1
        
        while time < max_steps:
            if random.random() < 0.6: # Probability p=0.6 to move to x+1
                state += 1
            else: # Probability q=0.4 to move to x-1
                state -= 1
            
            time += 1
            
            # Check for return to the set A = {0}
            if state == 0:
                returned_count += 1
                total_return_time += time
                break
    
    print("\n--- Simulation Results ---")
    print(f"Number of simulations started from state 0: {num_simulations}")
    print(f"Max steps before considering a path 'escaped': {max_steps}")
    print(f"Number of paths that returned to 0: {returned_count}")
    
    escaped_count = num_simulations - returned_count
    escape_percentage = (escaped_count / num_simulations) * 100
    
    print(f"Number of paths that escaped to infinity: {escaped_count} ({escape_percentage:.2f}%)")
    
    if returned_count > 0:
        avg_return_time = total_return_time / returned_count
        print(f"Average return time for paths that did return: {avg_return_time:.2f}")
    
    print("\n--- Conclusion ---")
    print("A significant fraction of paths did not return to the starting set within the maximum steps.")
    print("This suggests that the probability of returning is less than 1, meaning the chain is transient.")
    print("A transient chain cannot be positive recurrent. This aligns with our theoretical conclusion.")

simulate_non_positive_recurrent_chain()