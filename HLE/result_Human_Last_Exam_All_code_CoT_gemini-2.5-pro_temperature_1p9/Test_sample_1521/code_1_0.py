import random

def simulate_q_chain(start_pos, num_steps):
    """
    Simulates the h-transformed Markov chain (q-chain).

    The original p-chain is a symmetric random walk on {0, 1, 2, ...}, which is recurrent.
    A = {0}, h(x) = x.
    The resulting q-chain lives on {1, 2, ...}.
    Transitions are:
    - from i=1: must go to 2.
    - from i>1: go to i+1 with prob (i+1)/(2i), go to i-1 with prob (i-1)/(2i).
    """
    
    x = start_pos
    print(f"Starting simulation of the q-chain at position: {x}")
    print("--------------------------------------------------")
    
    for step in range(1, num_steps + 1):
        if x == 1:
            # From state 1, must move to state 2
            x = 2
        else:
            # For x > 1, calculate transition probabilities
            prob_up = (x + 1.0) / (2.0 * x)
            # prob_down = (x - 1.0) / (2.0 * x) is then 1 - prob_up
            
            if random.random() < prob_up:
                x = x + 1 # Move up
            else:
                x = x - 1 # Move down

        if step % (num_steps // 10) == 0:
            print(f"Step {step:6d}: Position = {x}")

if __name__ == '__main__':
    # Simulation parameters
    # The start position must be > 0
    START_POSITION = 1
    # Number of simulation steps
    NUMBER_OF_STEPS = 100000

    simulate_q_chain(START_POSITION, NUMBER_OF_STEPS)
    # The output illustrates the chain moving to larger and larger values,
    # suggesting transience.
    
    # Final answer for the two questions in the prompt
    final_answer = "(r, ?)"
    # The format requires printing each character including parenthesis, comma, space
    print("\nFinal answer to the theoretical questions:")
    print(f"The first character is: {final_answer[0]}")
    print(f"The second character is: {final_answer[1]}")
    print(f"The third character is: {final_answer[2]}")
    print(f"The fourth character is: {final_answer[3]}")
    print(f"The fifth character is: {final_answer[4]}")

