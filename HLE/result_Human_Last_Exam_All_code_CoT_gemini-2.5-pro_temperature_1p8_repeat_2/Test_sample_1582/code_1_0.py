import random
import statistics

def simulate_and_analyze_returns(target_state, num_returns):
    """
    Simulates a reflected random walk and analyzes the return times to a target state.
    
    This chain satisfies the conditions:
    - State space is {0, 1, 2, ...}
    - f(x) = x
    - A = {0}
    - E[f(X_1)|X_0=x] - f(x) = 0 >= 0 for x > 0.
    
    This chain is null recurrent, so we expect the average return time to grow.
    """
    
    print(f"Simulating a chain known to be null recurrent (not positive recurrent).")
    print(f"Calculating {num_returns} return times to state {target_state}...\n")
    
    current_state = target_state
    return_times = []
    steps_since_last_return = 0

    while len(return_times) < num_returns:
        steps_since_last_return += 1
        
        # Define transition probabilities for the reflected random walk
        if current_state == 0:
            # At state 0, must move to 1
            new_state = 1
        else:
            # At state x > 0, move to x-1 or x+1 with equal probability
            if random.random() < 0.5:
                new_state = current_state + 1
            else:
                new_state = current_state - 1
        
        current_state = new_state
        
        if current_state == target_state:
            return_times.append(steps_since_last_return)
            steps_since_last_return = 0

    print("--- Analysis of Return Times ---")
    print("If the chain were positive recurrent, the average return time would converge to a finite number.")
    print("Since it is null recurrent, we expect the running average to grow.\n")

    # Calculate and print cumulative averages at different points
    checkpoints = [10, 50, 100, 500, 1000]
    for cp in checkpoints:
        if cp <= num_returns:
            # Calculate the mean of the first 'cp' return times
            avg_return_time = statistics.mean(return_times[:cp])
            print(f"Average of the first {cp:4} return times: {avg_return_time:.2f}")

if __name__ == '__main__':
    # We will simulate 1000 returns to state 1
    # Note: This can take a few seconds as return times can be very long.
    simulate_and_analyze_returns(target_state=1, num_returns=1000)