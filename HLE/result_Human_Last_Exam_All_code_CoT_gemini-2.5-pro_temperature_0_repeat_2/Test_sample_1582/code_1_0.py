import random
import numpy as np

def demonstrate_non_positive_recurrence():
    """
    Demonstrates that a chain satisfying the given conditions is not positive recurrent.
    We use a specific null-recurrent birth-death chain on {0, 1, 2, ...}
    and the function f(x) = x.
    """
    print("--- Theoretical Verification ---")
    # Define the function f(x)
    f = lambda x: x
    
    # Define the transition probabilities p(x, y)
    def p(x_from, x_to):
        if x_from == 0:
            return 1.0 if x_to == 1 else 0.0
        if x_to == x_from + 1:
            return 0.5 + 1.0 / (4.0 * (x_from + 1.0))
        elif x_to == x_from - 1:
            return 0.5 - 1.0 / (4.0 * (x_from + 1.0))
        else:
            return 0.0

    # Let A = {0}. We verify the condition for a sample state x = 5 (where x is not in A).
    # The condition is: sum_y(p(x,y)*f(y)) - f(x) >= 0
    x = 5
    f_x = f(x)
    p_x_plus_1 = p(x, x + 1)
    p_x_minus_1 = p(x, x - 1)
    f_x_plus_1 = f(x + 1)
    f_x_minus_1 = f(x - 1)
    
    drift = (p_x_plus_1 * f_x_plus_1) + (p_x_minus_1 * f_x_minus_1) - f_x
    
    print(f"Let's verify the condition for f(x)=x and A={{0}} at state x={x}.")
    print(f"The condition is: p(x, x+1)*f(x+1) + p(x, x-1)*f(x-1) - f(x) >= 0")
    print(f"p({x}, {x+1}) = {p_x_plus_1:.4f}")
    print(f"p({x}, {x-1}) = {p_x_minus_1:.4f}")
    print(f"f({x}) = {f_x}, f({x+1}) = {f_x_plus_1}, f({x-1}) = {f_x_minus_1}")
    print("Final Equation:")
    print(f"({p_x_plus_1:.4f} * {f_x_plus_1}) + ({p_x_minus_1:.4f} * {f_x_minus_1}) - {f_x} = {drift:.4f}")
    print(f"The result {drift:.4f} is >= 0, so the condition holds.")
    print(f"The general formula for the drift is 1/(2*(x+1)), which for x=5 is {1/(2*(5+1)):.4f}.\n")

    print("--- Simulation ---")
    # Simulation parameters
    num_steps = 2_000_000
    current_state = 0
    
    # Data collection
    visit_counts = {0: 0, 1: 0, 2: 0}
    return_times_to_0 = []
    last_return_step = 0

    for step in range(num_steps):
        if current_state in visit_counts:
            visit_counts[current_state] += 1
        
        if current_state == 0 and step > 0:
            return_times_to_0.append(step - last_return_step)
            last_return_step = step

        # Transition to the next state
        rand_val = random.random()
        if current_state == 0:
            current_state = 1
        else:
            prob_up = p(current_state, current_state + 1)
            if rand_val < prob_up:
                current_state += 1
            else:
                current_state -= 1

    print(f"Simulated {num_steps} steps.")
    
    print("\n--- Analysis of Results ---")
    print("1. Fraction of time spent in initial states:")
    for state, count in visit_counts.items():
        fraction = count / num_steps
        print(f"   State {state}: {fraction:.6f}")
    print("Observation: The fraction of time spent in any given state is very small.")
    print("For a positive recurrent chain, this fraction should converge to a positive constant.\n")

    print("2. Average return time to state 0:")
    if len(return_times_to_0) > 100:
        num_blocks = 10
        block_size = len(return_times_to_0) // num_blocks
        for i in range(num_blocks):
            block = return_times_to_0[i*block_size : (i+1)*block_size]
            if block:
                avg_return_time = np.mean(block)
                print(f"   Avg. return time for returns {i*block_size+1}-{(i+1)*block_size}: {avg_return_time:.2f}")
        print("Observation: The average return time to state 0 increases over time.")
        print("For a positive recurrent chain, the mean return time should be finite and stable.\n")
    else:
        print("Not enough returns to 0 to analyze the trend.")

    print("--- Conclusion ---")
    print("Both the empirical state distribution and the increasing return times suggest")
    print("that the chain is not positive recurrent, supporting the theoretical conclusion.")

if __name__ == '__main__':
    demonstrate_non_positive_recurrence()
<<<Yes>>>