import math

# --- Step 1: Define assumptions and parameters ---

# Scenario 1: Synthesize a single Toffoli gate on a d=3 code
logical_t_gates_per_toffoli = 7

# Scenario 2: Run Shor's algorithm for a 2048-bit number on a d=5 code
# Number of Toffoli gates for Shor-2048 (Gidney & Ekerå, 2019)
toffoli_gates_for_shor_2048 = 2.5e7
logical_t_gates_for_shor_2048 = toffoli_gates_for_shor_2048 * logical_t_gates_per_toffoli

# Magic state distillation parameters
# Using the 15-to-1 protocol
distillation_input_size = 15

# --- Step 2: Calculation for Scenario 1 ---
# For a short demonstration (Toffoli gate), one level of distillation is sufficient.
distillation_levels_1 = 1
physical_t_gates_per_logical_t_1 = distillation_input_size ** distillation_levels_1
total_gates_scenario_1 = physical_t_gates_per_logical_t_1 * logical_t_gates_per_toffoli

# --- Step 3: Calculation for Scenario 2 ---
# For a long algorithm (Shor-2048), two levels of distillation are required
# to achieve a sufficiently low logical error rate.
distillation_levels_2 = 2
physical_t_gates_per_logical_t_2 = distillation_input_size ** distillation_levels_2
total_gates_scenario_2 = physical_t_gates_per_logical_t_2 * logical_t_gates_for_shor_2048

# --- Step 4: Final Calculation ---
total_non_clifford_gates = total_gates_scenario_1 + total_gates_scenario_2

# Print the final equation with all numbers
num1 = int(total_gates_scenario_1)
num2 = int(total_gates_scenario_2)
total_sum = int(total_non_clifford_gates)

print("The final calculation sums the physical non-Clifford T-gates from both scenarios.")
print("Equation:")
print(f"{num1} (from Scenario 1) + {num2} (from Scenario 2) = {total_sum}")
print("\nFinal Answer:")
print(f"{total_sum}")
