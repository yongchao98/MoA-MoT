import math

# This script calculates the minimal average resources required to simulate
# the correlations of a singlet quantum state using a local hidden variable model.

# Step 1: Define the CHSH values for the different models.
# The classical limit for any local hidden variable (LHV) model.
S_LHV = 2.0
# The quantum mechanical limit (Tsirelson's bound), achieved by a singlet state. This is our simulation target.
S_QUANTUM = 2 * math.sqrt(2)
# The limit for a non-local model using a single PR-Box or 1 bit of communication.
S_NONLOCAL = 4.0

# Step 2: Formulate the simulation as a probabilistic mixture.
# We want to find the minimal probability 'p' of using the non-local resource
# such that the average CHSH value matches the quantum limit.
# The governing equation is: S_QUANTUM = p * S_NONLOCAL + (1 - p) * S_LHV
# Solving for p: p = (S_QUANTUM - S_LHV) / (S_NONLOCAL - S_LHV)

p = (S_QUANTUM - S_LHV) / (S_NONLOCAL - S_LHV)

# Step 3: Display the results for each type of resource.

print("To simulate the correlations of a singlet state (which can reach a CHSH value of 2√2), a local model (limited to a CHSH of 2) must be supplemented with a non-local resource.\n")

# --- Resource 1: PR-Box ---
print("--- Simulation using PR-Boxes ---")
print("A non-signaling PR-Box can produce correlations with a CHSH value of 4.")
print("The simulation requires mixing the PR-Box strategy with a local one.")
print("The equation to satisfy is:")
# We use f-strings to embed the calculated values directly into the output.
# '1-p' is the probability of using the local model.
prob_local = 1 - p
print(f"  S_QUANTUM = p * S_PR_Box + (1-p) * S_LHV")
print(f"  {S_QUANTUM:.5f} = {p:.5f} * {S_NONLOCAL} + {prob_local:.5f} * {S_LHV}")
print(f"\nThis means the minimal average number of PR-Boxes required per trial is 'p'.")
print(f"Average PR-Boxes = {p:.8f}\n")


# --- Resource 2: Communication ---
print("--- Simulation using Classical Communication ---")
print("It is known that 1 bit of classical communication is sufficient to perfectly simulate a PR-Box, thus also achieving a CHSH value of 4.")
print("The calculation remains identical, where 'p' now represents the probability of using the 1-bit communication protocol.")
print("The equation to satisfy is:")
print(f"  S_QUANTUM = p * S_Communication(1bit) + (1-p) * S_LHV")
print(f"  {S_QUANTUM:.5f} = {p:.5f} * {S_NONLOCAL} + {prob_local:.5f} * {S_LHV}")
print(f"\nThe average communication cost is p * (1 bit) + (1-p) * (0 bits).")
print(f"Minimal average communication = {p:.8f} bits\n")

print(f"The required minimal resource cost for both cases is sqrt(2)-1.")
print(f"The calculated numerical value is {p}")
