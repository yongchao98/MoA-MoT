import random

def simulate_chain(num_returns):
    """
    Simulates a null-recurrent birth-death chain and calculates average return times.

    The chain is defined on {0, 1, 2, ...} with:
    - p(0, 1) = 1
    - p(x, x+1) = 0.5 for x > 0
    - p(x, x-1) = 0.5 for x > 0

    This chain satisfies the conditions in the problem description but is not
    positive recurrent. We demonstrate this by showing that the average return
    time to state 0 appears to grow indefinitely.
    """
    state = 0
    time = 0
    return_times = []
    total_time_at_return = 0

    print("Starting simulation to observe average return times to state 0.")
    print("-" * 60)
    
    # We start at state 0, so the first step is to leave.
    state = 1
    time = 1
    
    while len(return_times) < num_returns:
        # Move the chain one step
        if random.random() < 0.5:
            state += 1
        else:
            state -= 1
        time += 1

        # Check for return to state 0
        if state == 0:
            return_times.append(time)
            total_time_at_return += time
            current_avg = total_time_at_return / len(return_times)
            
            # Print the running average at certain intervals
            if len(return_times) % (num_returns // 10 or 1) == 0:
                print(f"Return #{len(return_times):<5}: Current time = {time:<10} | Average return time = {current_avg:.2f}")

            # Reset for the next excursion from 0
            # To measure time between returns, we reset the timer
            # For total time, we don't. Let's use total time to see growth.
            # Next excursion starts now
            time = 0
            state = 1
            time = 1


if __name__ == '__main__':
    # We will simulate for a number of returns to state 0.
    # A larger number gives a better idea, but takes longer.
    # Note: For a null-recurrent chain, return times can be very long.
    simulate_chain(num_returns=100)
