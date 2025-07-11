import numpy as np

# --- Problem Setup and Theoretical Verification ---
# Let's define a Markov chain and a function f that satisfy the problem's conditions.

# 1. The Markov Chain: Biased Random Walk on {0, 1, 2, ...}
# Let the probability of moving right (x -> x+1) be p and left (x -> x-1) be q.
# We choose p > q to create a drift to infinity, making the chain not positive recurrent.
p = 0.6
q = 1 - p
# The transitions are p(x, x+1) = p, p(x, x-1) = q for x >= 1.
# At state 0, let's assume the particle must move to state 1, so p(0,1) = 1.
# This chain is irreducible and is known to be transient, thus NOT positive recurrent.

# 2. The Finite Set A and the Function f
A = {0}
# We choose the function f(x) = x.
# f is non-negative on the state space and f(x) -> infinity as x -> infinity.
def f(x):
    return x

# 3. Verification of the Drift Condition
# We check if E[f(X_{n+1}) | X_n = x] - f(x) >= 0 for all x not in A (i.e., x >= 1).
# The expected value of f at the next step is p*f(x+1) + q*f(x-1).
# The drift is [p*f(x+1) + q*f(x-1)] - f(x)
# which for f(x)=x is [p*(x+1) + q*(x-1)] - x = px + p + qx - q - x = (p+q)x + p - q - x = p - q.
drift = p - q

print("--- Verification of Conditions ---")
print(f"Chain: Biased random walk with p(x,x+1)={p} for x>=1.")
print(f"Set A = {A}, Function f(x) = x.")
print("The drift E[f(X_n+1) | X_n=x] - f(x) for x not in A is constant.")
print("Equation for drift: p * f(x+1) + q * f(x-1) - f(x)")
print("Substituting f(x)=x: p * (x+1) + q * (x-1) - x")
print("This simplifies to: p - q")
print(f"For our values, the drift calculation is:")
print(f"{p} - {q} = {drift:.2f}")
print(f"Since {drift:.2f} >= 0, the condition is satisfied.")


# --- Numerical Simulation to Illustrate the Result ---
print("\n--- Numerical Simulation ---")
num_walks = 5000
max_steps = 200
start_pos = 10  # Starting outside A

print(f"Simulating {num_walks} walks for {max_steps} steps, starting at X_0 = {start_pos}.")
print("We will track the average position, which is the average of f(X_n).")

# Initialize all walks at the starting position
positions = np.full(num_walks, start_pos, dtype=int)
report_steps = [10, 50, 100, 150, 200]

for step in range(1, max_steps + 1):
    # Determine moves for walks at positions > 0
    at_positive = positions > 0
    num_positive = np.sum(at_positive)
    
    # Generate random moves [1, -1] with probabilities [p, q]
    moves = np.zeros(num_walks, dtype=int)
    if num_positive > 0:
      moves[at_positive] = np.random.choice([1, -1], size=num_positive, p=[p, q])

    # Walks that are at state 0 must move to 1
    at_zero = positions == 0
    moves[at_zero] = 1

    # Update all positions
    positions += moves

    # Report the average position at specified steps
    if step in report_steps:
        # The average of f(X_n) is the average position
        avg_f_value = np.mean(positions)
        print(f"Step {step:3d}: Average f(X_n) = Average Position = {avg_f_value:.2f}")

print("\nSimulation Conclusion: The average value of f(X_n) consistently increases,")
print("illustrating the chain's tendency to drift away from A towards infinity.")
