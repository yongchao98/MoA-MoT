import random
import math

def simulate_escape(n, num_simulations=10000):
    """
    Simulates a simple random walk to estimate the escape probability p_n.
    A walk 'escapes' if it leaves the cube C_n before returning to the start.
    """
    escapes = 0
    start_pos = (n, 0, 0)
    
    for _ in range(num_simulations):
        pos = start_pos
        
        # First step
        moves = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
        move = random.choice(moves)
        pos = (pos[0] + move[0], pos[1] + move[1], pos[2] + move[2])

        while True:
            # Check for return to start
            if pos == start_pos:
                break # Failed to escape
            
            # Check for escape from C_n
            if not (0 <= pos[0] <= 2*n and 0 <= pos[1] <= 2*n and 0 <= pos[2] <= 2*n):
                escapes += 1
                break
                
            # Take another step
            move = random.choice(moves)
            pos = (pos[0] + move[0], pos[1] + move[1], pos[2] + move[2])
            
    return escapes / num_simulations

def main():
    """
    Calculates the limit by simulating for increasing n and printing the results.
    """
    # Using a range of n values to observe the trend
    n_values = [5, 10, 20, 40, 80]
    
    print("This program numerically verifies the limit value.")
    print("It simulates the random walk for several values of n.")
    print("The final column should approach 1 as n increases.")
    print("-" * 55)
    print(f"{'n':>5s} {'p_n (estimated)':>20s} {'ln(1/p_n)/ln(n)':>25s}")
    print("-" * 55)

    # Note: num_simulations should be large for good accuracy, but will take longer.
    # For a quick run, 10,000 is okay. For better estimates, use 100,000 or more.
    num_simulations_per_n = 20000 

    for n in n_values:
        p_n = simulate_escape(n, num_simulations=num_simulations_per_n)
        
        if p_n > 0:
            limit_val = math.log(1 / p_n) / math.log(n)
            print(f"{n:>5d} {p_n:>20.6f} {limit_val:>25.6f}")
        else:
            print(f"{n:>5d} {'p_n is too small to measure':>20s} {'N/A':>25s}")
    
    print("-" * 55)
    print("The theoretical value of the limit is 1.")

if __name__ == "__main__":
    main()
