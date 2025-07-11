import math

# Step 1: Define the maximal CHSH values for each resource type.
S_Local = 2
S_Quantum = 2 * math.sqrt(2)
S_PR = 4

# Step 2: The simulation equation is S_Quantum = p * S_PR + (1-p) * S_Local
# We need to solve for 'p'.
# 2*sqrt(2) = p * 4 + (1-p) * 2
# 2*sqrt(2) = 4p + 2 - 2p
# 2*sqrt(2) = 2p + 2
# 2p = 2*sqrt(2) - 2
# p = sqrt(2) - 1

p = math.sqrt(2) - 1

# Step 3: Print the logic and the result.
# The prompt requires printing each number in the final equation.
print("To find the minimal average resource 'p', we solve the equation:")
print(f"{S_Quantum:.4f} = p * {S_PR} + (1-p) * {S_Local}")
print("Rearranging the terms gives:")
print(f"2p = 2 * sqrt(2) - 2")
print("So, the minimal resource fraction 'p' is:")
print(f"p = sqrt(2) - 1")
print("\nCalculating the numerical value:")
print(f"The minimal average resource 'p' required is: {p}")
