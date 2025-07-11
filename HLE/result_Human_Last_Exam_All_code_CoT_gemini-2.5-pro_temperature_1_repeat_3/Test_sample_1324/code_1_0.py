import math

def simulate_chocolate_game(initial_chocolates, max_steps):
    """
    Simulates the chocolate passing game and prints the state at each step.

    Args:
        initial_chocolates (list): A list of even integers representing the initial
                                   number of chocolates for each person.
        max_steps (int): The number of minutes to simulate.
    """
    chocolates = list(initial_chocolates)
    n = len(chocolates)

    print(f"Starting simulation with n={n} people.")
    print(f"Initial state C^0: {chocolates}\n")

    for i in range(max_steps):
        # Calculate max, min, and difference for the current step
        h_i = max(chocolates)
        l_i = min(chocolates)
        d_i = h_i - l_i

        print(f"--- Minute i = {i} ---")
        print(f"Chocolates C^{i}: {chocolates}")
        print(f"Max h^{i} = {h_i}, Min l^{i} = {l_i}, Diff d^{i} = {d_i}")

        # If the system has stabilized, we can stop
        if d_i == 0:
            print("\nSystem has reached equilibrium. All chocolate counts are equal.")
            print("For this step (and all future steps):")
            print(f"Statement 1 fails: We need d^{i+m} < d^{i}, but d^{i+m}=0 and d^{i}=0. (0 < 0 is false)")
            print(f"Statement 3 fails: We need l^{i+m} > l^{i}, but l^{i+m}={l_i} and l^{i}={l_i}. ({l_i} > {l_i} is false)")
            # We can also check statement 2 by observing the sequence of l_i values
            # which will be non-decreasing.
            break

        # Calculate the next state
        next_chocolates = [0] * n
        for k in range(n):
            # Person k receives from person k-1 (with wrap-around)
            person_k_chocs = chocolates[k]
            neighbor_chocs = chocolates[(k - 1 + n) % n]

            # Each person keeps half and receives half from the neighbor
            new_total = (person_k_chocs / 2) + (neighbor_chocs / 2)

            # If the result is odd, take an extra one
            if new_total % 2 != 0:
                new_total += 1

            next_chocolates[k] = int(new_total)

        chocolates = next_chocolates
        print("")


# Example from the thought process: n=4, initial chocolates [10, 2, 2, 2]
# This setup clearly shows the dynamics.
# All initial values must be even.
initial_setup = [10, 2, 2, 2]
simulation_steps = 10
simulate_chocolate_game(initial_setup, simulation_steps)
