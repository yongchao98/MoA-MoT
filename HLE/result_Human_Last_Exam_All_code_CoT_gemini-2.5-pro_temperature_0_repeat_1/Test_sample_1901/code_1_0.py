import math

# 1. Define the target and initial parameters based on the plan
T_target = 22.5  # °C, target room temperature transition
n_initial = 5      # Initial chain length (pentyl)
T_initial = 35.0   # °C, known transition temp for 4-pentylcyanobenzene (n=5)
k = 4.0        # °C, assumed temperature decrease per added CH2 group

# 2. Solve for the final chain length 'n'
# Equation: T_target = T_initial - k * (n_final - n_initial)
# Rearranged: n_final = n_initial + (T_initial - T_target) / k
n_final_float = n_initial + (T_initial - T_target) / k

# Round to the nearest integer for the final chain length
n_final_int = int(round(n_final_float))

# 3. Output the process and the final equation with numbers
print("Objective: Find the alkyl chain length 'n' for a CnH2n+1-Ph-CN molecule with a transition temperature near room temperature.")
print(f"Target Temperature (T_target): {T_target}°C")
print("-" * 30)
print("Modeling based on the provided rules:")
print(f"Initial molecule has chain length n_initial = {n_initial}.")
print(f"Its known transition temperature is T_initial = {T_initial}°C.")
print(f"We assume the temperature decreases by k = {k}°C for each carbon added to the chain.")
print("-" * 30)
print("We use the equation: T_target = T_initial - k * (n_final - n_initial)")
print("To find the required chain length 'n_final', we solve for it.")
print("\nHere is the final equation with the numbers plugged in:")
# Output each number in the final equation as requested.
print(f"{T_target} = {T_initial} - {k} * (n_final - {n_initial})")

print(f"\nSolving for n_final gives a value of {n_final_float:.2f}.")
print(f"Since the chain length 'n' must be an integer, we round this to the nearest whole number.")
print(f"\nThe recommended integer chain length is n = {n_final_int}.")
print(f"The suggested molecule is 4-octylcyanobenzene.")