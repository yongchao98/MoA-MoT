import sympy

# Define symbols
v_L, p_L, v_H, p_H = sympy.symbols('v_L p_L v_H p_H')

# Based on the derivation, we found the equilibrium contracts.
# Low-type contract
eq_v_L = 1
eq_p_L = 1

# High-type contract
eq_v_H = 30
eq_p_H = 30

# Now, we print the final answer and verify the equilibrium conditions.

print("Separating Equilibrium Contracts:")
print(f"Low-Type Contract (v_L, p_L): ({eq_v_L}, {eq_p_L})")
print(f"High-Type Contract (v_H, p_H): ({eq_v_H}, {eq_p_H})")
print("\n--- Verification of Equilibrium Conditions ---")

# 1. Productivity Constraints
v_L_max = 15 # Solved from v_L <= 20 - (1/3)*v_L
v_H_max = 30
print("\n1. Productivity Constraints:")
print(f"v_L <= {v_L_max}? {eq_v_L <= v_L_max}")
print(f"v_H <= {v_H_max}? {eq_v_H <= v_H_max}")


# 2. Individual Rationality (IR)
u_L_contract = eq_p_L - eq_v_L**2
u_H_contract = eq_p_H - eq_v_H
print("\n2. Individual Rationality (IR):")
print(f"Low type utility >= 0? {u_L_contract} >= 0 -> {u_L_contract >= 0}")
print(f"High type utility >= 0? {u_H_contract} >= 0 -> {u_H_contract >= 0}")

# 3. Incentive Compatibility (IC)
# Low type's utility from choosing High type's contract
u_L_mimic = eq_p_H - eq_v_H**2
# High type's utility from choosing Low type's contract
u_H_mimic = eq_p_L - eq_v_L

print("\n3. Incentive Compatibility (IC):")
print("Does Low type prefer (v_L, p_L)?")
print(f"  Utility from own contract: {u_L_contract}")
print(f"  Utility from mimicking High type: {u_L_mimic}")
print(f"  Is {u_L_contract} >= {u_L_mimic}? {u_L_contract >= u_L_mimic}")

print("Does High type prefer (v_H, p_H)?")
print(f"  Utility from own contract: {u_H_contract}")
print(f"  Utility from mimicking Low type: {u_H_mimic}")
print(f"  Is {u_H_contract} >= {u_H_mimic}? {u_H_contract >= u_H_mimic}")
