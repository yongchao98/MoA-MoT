import math

# Introduction to the problem and resources
print("To simulate the quantum correlations from any POVM measurements on a singlet state")
print("using a Local Hidden Variable (LHV) model, supplemental resources are required.")
print("The two primary resources considered are classical communication and non-signaling PR-boxes.")
print("Here, we calculate the minimal average amount required for each type of resource.")
print("-" * 65)

# --- Resource 1: Classical Communication ---
print("Resource 1: Minimal Average Classical Communication (C)")
print("The Toner-Bacon protocol demonstrates that the correlations of a singlet state")
print("can be perfectly simulated with an LHV model supplemented by classical communication.")
communication_cost = 1
print(f"The proven minimal average cost is: C = {communication_cost} bit")
print("-" * 65)

# --- Resource 2: Non-Signaling PR-Box ---
print("Resource 2: Minimal Average PR-Box Fraction (q)")
print("This cost is the minimal fraction 'q' of a PR-box needed in a mixture")
print("with a purely local resource to reproduce the quantum correlations.")
print("\nThe simulation must match the quantum maximum for the CHSH inequality,")
print("known as Tsirelson's bound.")
print("\nThe governing equation is: q * S_PR + (1 - q) * S_L = S_QM")

S_PR = 4
S_L = 2
val_sqrt2 = math.sqrt(2)
S_QM_val = 2 * val_sqrt2

print(f"\nWhere:")
print(f"S_PR (PR-Box CHSH Score) = {S_PR}")
print(f"S_L  (Local CHSH Limit)   = {S_L}")
print(f"S_QM (Quantum CHSH Limit) = 2 * sqrt(2) \u2248 {S_QM_val:.6f}")

print("\nSolving the equation for q:")
print(f"q * {S_PR} + (1 - q) * {S_L} = 2 * sqrt(2)")
print(f"{S_PR}q + {S_L} - {S_L}q = 2 * sqrt(2)")
print(f"{(S_PR - S_L)}q + {S_L} = 2 * sqrt(2)")
print(f"2q + 2 = 2 * sqrt(2)")
print( "q + 1 = sqrt(2)")
print( "q = sqrt(2) - 1")

# Final calculation
pr_box_fraction = val_sqrt2 - 1

print("\nFinal Calculation:")
print(f"q = {val_sqrt2} - 1")
print(f"q \u2248 {pr_box_fraction}")
print("-" * 65)

print("\nSummary of Minimal Average Resources:")
print(f"1. Classical Communication Cost: {communication_cost} bit")
print(f"2. PR-Box Resource Cost:       sqrt(2) - 1 \u2248 {pr_box_fraction:.6f}")
print("-" * 65)
<<<0.4142135623730951>>>