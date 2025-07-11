import random

def simulate_transient_chain(steps, initial_state=0):
    """
    Simulates a Markov chain on {0, 1, 2, ...} that satisfies the conditions
    and demonstrates it is not positive recurrent.

    The chain is defined by the following transition probabilities:
    - p(0, 1) = 1
    - For x > 0, p(x, x+1) = 0.6 and p(x, x-1) = 0.4

    This chain is irreducible on the state space {0, 1, 2, ...}.
    Let's check the conditions from the problem:
    - Let the function f be f(x) = x. f is non-negative and f(x) -> infinity.
    - Let the finite set A be A = {0}.
    - For any state x not in A (i.e., x > 0), we check the condition:
      E[f(X_{n+1}) | X_n=x] - f(x)
      = E[X_{n+1} | X_n=x] - x
      = (0.6 * (x+1) + 0.4 * (x-1)) - x
      = (0.6*x + 0.6 + 0.4*x - 0.4) - x
      = (x + 0.2) - x
      = 0.2
    Since 0.2 >= 0, the condition is met.
    
    Based on the theory, this chain should not be positive recurrent. The simulation
    will show the chain drifting to infinity.
    """
    current_state = initial_state
    path = [current_state]

    print(f"Simulating a Markov chain for {steps} steps, starting at state {initial_state}.")
    print("This chain satisfies the conditions and is expected to be transient (drifting to +infinity).")
    print("-" * 40)

    # First, run the simulation and print the state at intervals to show the drift
    for i in range(steps):
        if current_state == 0:
            current_state = 1
        else:
            # random() returns a float in [0.0, 1.0)
            if random.random() < 0.6:
                current_state += 1
            else:
                current_state -= 1
        
        # Print the state at regular intervals to observe the drift
        if (i + 1) % (steps // 10) == 0:
            print(f"Step {i+1:7d}: Current state = {current_state}")
            
    print("-" * 40)
    print(f"Final state after {steps} steps: {current_state}")
    
    # In a positive recurrent chain, the chain returns to the origin frequently.
    # To quantify this, we re-run the simulation and count visits to the origin.
    current_state = initial_state
    visits_to_origin = 1 if initial_state == 0 else 0
    for i in range(steps):
        if current_state == 0:
            current_state = 1
        else:
            if random.random() < 0.6:
                current_state += 1
            else:
                current_state -= 1
        if current_state == 0:
            visits_to_origin += 1
            
    print(f"Number of visits to the origin in {steps} steps: {visits_to_origin}")
    if visits_to_origin > (1 if initial_state == 0 else 0):
        print("The chain returned to the origin.")
    else:
        print("The chain did not return to the origin after leaving it.")
    print("\nAs the number of steps increases, the chain tends to move further away,")
    print("making returns to the origin increasingly rare. This demonstrates behavior")
    print("that is not positive recurrent.")

# Execute the simulation for a large number of steps.
simulate_transient_chain(steps=100000)