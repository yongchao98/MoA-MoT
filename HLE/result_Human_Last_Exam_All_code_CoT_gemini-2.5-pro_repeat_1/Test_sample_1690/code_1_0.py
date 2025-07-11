import math

# Step 1: Define initial parameters
P = 100  # Initial power in W
x = 4    # Control subsystem power in W
y = 5    # Camera subsystem power in W
z = 1    # Sensor subsystem power in W
T_half = 400  # Half-life in days

# Step 2: Calculate shutdown times for a and b
# The formula to find the time t is: t = T_half * (ln(Target_Power) - ln(P)) / ln(0.5)

# a. Time to turn off the camera subsystem
target_power_a = x + y + z
t_a = T_half * (math.log(target_power_a) - math.log(P)) / math.log(0.5)
days_a = round(t_a)

# b. Time to turn off the sensor subsystem
target_power_b = x + z
t_b = T_half * (math.log(target_power_b) - math.log(P)) / math.log(0.5)
days_b = round(t_b)

# Step 3: Calculate memory usage for c
# Data type sizes: int = 5D, frac = 6D
# A memory-efficient program would declare variables for:
# - Inputs: P, x, y, z (4 ints)
# - Results: days_a, days_b (2 ints)
# - Time-saving pre-calculations: ln_P, ln_half (2 fracs)
# - Temporary calculation variables: temp_ln, result (2 fracs)
# Total: 6 int variables, 4 frac variables
size_int_D = 5
size_frac_D = 6
num_ints = 6
num_fracs = 4
memory_c = (num_ints * size_int_D) + (num_fracs * size_frac_D)

# Step 4: Calculate the number of exp function calls for d
# The program must calculate ln(0.5), ln(P), ln(x+y+z), and ln(x+z).
# Number of exp calls for ln(v) is ceil(log2(v * 1e9)).
values_for_ln = [0.5, P, target_power_a, target_power_b]
exp_calls_d = 0
for val in values_for_ln:
    # k >= log2(val * 1e9)
    k = math.log2(val * 1e9)
    exp_calls_d += math.ceil(k)

# Final output in the required format
print(f"{days_a}:{days_b}:{memory_c}:{exp_calls_d}")