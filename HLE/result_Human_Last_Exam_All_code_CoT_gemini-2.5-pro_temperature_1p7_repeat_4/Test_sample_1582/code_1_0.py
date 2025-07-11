import random

def simulate_srw(num_steps):
    """
    Simulates a simple symmetric random walk on N_0 with a reflective barrier at 0.
    
    Args:
        num_steps (int): The total number of steps to simulate.
        
    Returns:
        float: The proportion of time the chain spent at state 0.
    """
    state = 0
    visits_to_zero = 0
    
    for _ in range(num_steps):
        if state == 0:
            visits_to_zero += 1
            # Transition from 0 is always to 1
            state = 1
        else:
            # Transition from x > 0 is to x-1 or x+1 with probability 0.5
            if random.random() < 0.5:
                state += 1
            else:
                state -= 1
    
    return visits_to_zero / num_steps

def main():
    """
    Main function to run the simulation and demonstrate the property.
    """
    print("This script simulates a Markov chain that satisfies the conditions in the problem.")
    print("The chain is a simple symmetric random walk on {0, 1, 2, ...} with reflection at 0.")
    print("The conditions are satisfied with A={0} and f(x)=x.")
    print("-" * 50)

    # Let's check the condition Pf(x) >= f(x) for a sample x, with f(x)=x
    x = 5
    f_x = x
    p_x_plus_1 = 0.5
    f_x_plus_1 = x + 1
    p_x_minus_1 = 0.5
    f_x_minus_1 = x - 1
    
    pf_x = p_x_plus_1 * f_x_plus_1 + p_x_minus_1 * f_x_minus_1
    drift = pf_x - f_x

    print(f"Checking the drift condition for x = {x} with f(x)=x:")
    print(f"E[f(X_1)|X_0={x}] - f({x}) = ({p_x_plus_1} * {f_x_plus_1} + {p_x_minus_1} * {f_x_minus_1}) - {f_x} = {pf_x} - {f_x} = {drift}")
    print(f"Since the result {drift} is >= 0, the condition holds.")
    print("-" * 50)

    print("Running simulations to find the proportion of time spent at state 0.")
    print("For a chain that is not positive recurrent, this proportion should tend to 0.")
    
    # We use a set of increasing simulation steps to observe the trend
    step_counts = [10000, 50000, 100000, 500000, 1000000, 5000000]
    
    for steps in step_counts:
        # For reproducibility of the output
        random.seed(0)
        proportion = simulate_srw(steps)
        print(f"For {steps: >7} steps, the proportion of time at state 0 is: {proportion:.6f}")
        
    print("\nAs the number of steps increases, the proportion of time at state 0 decreases.")
    print("This is the expected behavior for a chain that is not positive recurrent.")

if __name__ == "__main__":
    main()