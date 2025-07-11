import random

def simulate_chain():
    """
    Simulates a Markov chain on Z that satisfies the given conditions
    and is not positive recurrent.
    """
    x = 0  # Start at the origin
    total_steps = 20000
    
    print("Simulating a Markov chain that is not positive recurrent.")
    print(f"Starting position X_0: {x}")
    print("-" * 30)

    for n in range(1, total_steps + 1):
        # Determine the probability of moving to x+1
        if x > 10:
            # For x > 10, there's a drift pushing it further positive.
            # p(x, x+1) = 0.6
            p_plus_1 = 0.6
        elif x < -10:
            # For x < -10, there's a drift pushing it further negative.
            # p(x, x-1) = 0.6, so p(x, x+1) = 0.4
            p_plus_1 = 0.4
        else: # x is in [-10, 10]
            # Inside the set A, the walk is a standard symmetric random walk.
            p_plus_1 = 0.5
        
        # Perform one step of the Markov chain
        if random.random() < p_plus_1:
            x += 1
        else:
            x -= 1
        
        # Print the position at certain intervals to observe its behavior
        if n % 2000 == 0:
            print(f"Position at step n={n}: X_n = {x}")

    print("-" * 30)
    print("Notice how the position tends to drift away from the origin.")
    print("This illustrates the chain is not positive recurrent; it does not return")
    print("to the origin 'often enough'.")

if __name__ == "__main__":
    simulate_chain()