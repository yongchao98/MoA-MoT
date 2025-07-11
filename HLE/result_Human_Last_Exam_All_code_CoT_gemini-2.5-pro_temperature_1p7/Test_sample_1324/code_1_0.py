import math

def simulate_chocolate_passing(initial_chocolates, max_steps=20):
    """
    Simulates the chocolate passing game.

    Args:
        initial_chocolates (list): A list of even integers representing initial chocolates.
        max_steps (int): The maximum number of simulation steps.
    """
    n = len(initial_chocolates)
    chocolates = list(initial_chocolates)

    print(f"Starting simulation with n={n} people.")
    print(f"Initial state c^0: {initial_chocolates}\n")

    for i in range(max_steps):
        # Calculate and print stats for the current state
        h_i = max(chocolates)
        l_i = min(chocolates)
        d_i = h_i - l_i
        
        print(f"Minute i={i}:")
        print(f"  Chocolates: {chocolates}")
        print(f"  h^{i}={h_i}, l^{i}={l_i}, d^{i}={d_i}")
        
        # If all have the same number, the process stops
        if d_i == 0:
            print("\nAll people have the same number of chocolates. The system has stabilized.")
            final_equation = " = ".join(map(str, chocolates))
            print(f"Final equation: {final_equation}")
            return

        # Prepare for the next minute
        next_chocolates = [0] * n
        chocolates_added = 0
        for k in range(n):
            # Get the neighbor who passes chocolates to person k
            # In the problem, p_k gives to p_{k+1}, so p_k receives from p_{k-1}.
            # The problem defines c_b^{i+1} based on c_b^i and c_a^i, where a passes to b.
            # So for person k, a=k-1 and b=k.
            # (Indices are 0 to n-1, so we use modulo arithmetic)
            c_self = chocolates[k]
            c_neighbor = chocolates[(k - 1 + n) % n]
            
            # Each person passes half, so person k receives c_neighbor/2 and keeps c_self/2
            total = (c_self + c_neighbor) / 2
            
            # If odd, they get an extra one
            if total % 2 != 0:
                next_chocolates[k] = math.ceil(total)
                chocolates_added +=1
            else:
                next_chocolates[k] = int(total)

        chocolates = next_chocolates
        print(f"  (Chocolates added this round: {chocolates_added})\n")
        
    print("\nMax steps reached. The system may not have stabilized.")
    print(f"Final state: {chocolates}")

if __name__ == '__main__':
    # Example from the problem description's logic derivation
    # n=5, two blocks of people
    initial_config = [10, 2, 2, 2, 10]
    simulate_chocolate_passing(initial_config)
