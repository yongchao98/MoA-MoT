import random

def simulate_and_analyze():
    """
    Simulates a Markov chain that satisfies the problem's conditions and
    demonstrates that it is not positive recurrent.

    The example chain is a random walk on non-negative integers (Z_+):
    - State space Sigma = {0, 1, 2, ...}
    - From x > 0, move to x+1 with probability p, and to x-1 with probability 1-p.
    - From x = 0, move to 1 with probability 1 (reflecting barrier).

    Let f(x) = x and the finite set be A = {0}.
    - f(x) is non-negative and f(x) -> infinity as x -> infinity.
    - We check the condition: E[f(X_n+1)|X_n=x] - f(x) >= 0 for x not in A (x > 0).
      The calculation is: p*(x+1) + (1-p)*(x-1) - x = 2p - 1.
      This is non-negative if p >= 0.5.

    We choose p = 0.6, so the condition holds. This chain is transient. We will
    show it's not positive recurrent by simulating it and observing that the
    fraction of visits to state 0 approaches zero over time.
    """

    # Parameters for the simulation
    p = 0.6  # Probability of moving right (p > 0.5)
    num_steps = 2000000

    print(f"Simulating a Markov chain with p = {p} for {num_steps} steps.")
    print("This setup satisfies the conditions given in the problem.")
    print("If the chain is not positive recurrent, the fraction of visits to any state (e.g., 0) should decay to 0.")
    print("-" * 60)

    # Simulation state
    current_state = 0
    visit_count_state_0 = 1  # We start at state 0

    # Reporting setup
    report_points = [int(num_steps * i / 10) for i in range(1, 11)]

    print(f"{'Total Steps':>15} | {'Visits to State 0':>20} | {'Fraction at State 0':>20}")
    print("-" * 60)

    for step in range(1, num_steps + 1):
        # Move to the next state
        if current_state == 0:
            current_state = 1
        else:
            if random.random() < p:
                current_state += 1
            else:
                current_state -= 1
        
        # Count visits
        if current_state == 0:
            visit_count_state_0 += 1

        # Print report at intervals
        if step in report_points:
            fraction = visit_count_state_0 / step
            print(f"{step:>15,} | {visit_count_state_0:>20,} | {fraction:>20.6f}")

    print("-" * 60)
    print("As observed, the fraction of time spent at state 0 decreases as the number of steps increases,")
    print("which is characteristic of a chain that is not positive recurrent.\n")

    # Outputting the numbers in the final equation as requested.
    print("Final check of the condition E[f(X_n+1)|X_n=x] - f(x) >= 0:")
    # The left side of the equation can be calculated for our f(x)=x and p=0.6
    p_val = 0.6
    q_val = 1 - p_val
    drift = 2 * p_val - 1
    
    print(f"f(x) = x, A = {{0}}, p = {p_val}")
    print(f"For x > 0, E[f(X_n+1)|X_n=x] - f(x) = {p_val}*(x+1) + {q_val:.1f}*(x-1) - x")
    print(f"                                   = ({p_val} + {q_val:.1f})*x + {p_val} - {q_val:.1f} - x")
    print(f"                                   = 1*x + ({p_val - q_val:.1f}) - x")
    print(f"                                   = {drift:.1f}")
    print(f"The result {drift:.1f} is >= 0, so the condition holds.")

# Run the simulation
simulate_and_analyze()
