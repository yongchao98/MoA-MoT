import math

# The problem is to find the minimal resources (PR-Boxes, Communication) needed for a
# Local Hidden Variable (LHV) model to simulate the correlations of a singlet quantum state
# for any set of POVM measurements.

# According to the resource theory of non-locality, simulating the correlations of one
# maximally entangled state (a singlet state, or 'ebit') requires a combination of
# a non-local resource (a PR-Box) and classical communication.

# The amount of one-way communication 'C' required can be calculated. For a pure
# entangled state |ψ⟩ = cos(θ)|00⟩ + sin(θ)|11⟩, the cost is C = -log₂(cos²(θ)).
# A singlet state is maximally entangled, which corresponds to θ = π/4.

# We calculate the communication cost C. The number of PR-boxes is known to be 1 from the theory.
num_pr_box = 1.0
num_singlet_state = 1.0

# Define the angle for a maximally entangled state in radians.
theta = math.pi / 4

# Calculate the communication cost using the formula.
# communication_cost = -math.log2( (math.cos(math.pi/4))**2 )
# = -math.log2( (1/math.sqrt(2))**2 )
# = -math.log2(1/2)
# = 1.0
communication_cost = -math.log2(math.cos(theta)**2)


print("The minimal resources required to simulate a singlet state with an LHV model are given by the equivalence relation:")
print("-" * 70)
print(f"{num_singlet_state} * [Singlet State]  =  {num_pr_box} * [PR-Box]  +  {communication_cost} * [bit of one-way communication]")
print("-" * 70)
print("\nThis means that on average, for each pair of measurements, the simulation requires:")
print(f"- {num_pr_box} use of a non-signaling PR-Box.")
print(f"- {communication_cost} bit of one-way classical communication.")
