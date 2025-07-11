import math

# --- Step 1: Define the known values and constants ---

# Activity at the first measurement (kBq/mL)
A1 = 1.4
# Activity at the second measurement (kBq/mL)
A2 = 2.1
# Time between the two measurements (days)
delta_t_days = 14

# Half-life of Yttrium-90 (hours)
t_half_Y_hours = 64.1

# --- Step 2: Calculate the decay constant for Y-90 in units of days^-1 ---

# Convert half-life from hours to days
t_half_Y_days = t_half_Y_hours / 24
# Calculate the decay constant (lambda = ln(2) / half_life)
lambda_Y_per_day = math.log(2) / t_half_Y_days

# --- Step 3: Set up and solve the equations ---

# We have the system of equations:
# 1) A1 = A_Sr * (1 - exp(-lambda * T))
# 2) A2 = A_Sr * (1 - exp(-lambda * (T + delta_t_days)))
#
# Dividing equation 1 by equation 2 gives:
# A1 / A2 = [1 - exp(-lambda * T)] / [1 - exp(-lambda * T) * exp(-lambda * delta_t_days)]
#
# Let's define x = exp(-lambda * T) and k = exp(-lambda * delta_t_days)
# A1 / A2 = (1 - x) / (1 - x*k)
#
# Rearranging to solve for x:
# x = 1 / ((A2/A1)*(1-k) + k) --> A less intuitive form
# Let's rearrange as in the plan:
# (A1/A2) * (1 - x*k) = 1 - x
# A1/A2 - (A1/A2)*x*k = 1 - x
# x - (A1/A2)*x*k = 1 - (A1/A2)
# x * (1 - (A1/A2)*k) = 1 - A1/A2
# x = (1 - A1/A2) / (1 - (A1/A2)*k) --> This is correct and stable.

# Calculate the constant k = exp(-lambda * delta_t_days)
k = math.exp(-lambda_Y_per_day * delta_t_days)

# Calculate the ratio of activities
ratio = A1 / A2

# Solve for x
# ratio = (1-x) / (1-x*k) => ratio*(1-x*k) = 1-x => ratio - ratio*x*k = 1 - x => x(1 - ratio*k) = 1-ratio => x = (1-ratio)/(1-ratio*k)
x = (1 - ratio) / (1 - ratio * k)

# Now, solve for T using x = exp(-lambda * T)
# -lambda * T = ln(x)
# T = -ln(x) / lambda
T_days = -math.log(x) / lambda_Y_per_day

# --- Step 4: Print the final answer and the equation ---

print("This script solves for the time (T) between sample irradiation and the first measurement.")
print("The governing equation is derived from the Y-90 ingrowth formula:")
print(f"A(t) = A_Sr * (1 - exp(-lambda_Y * t))")
print("\nWe have a system of two equations:")
print(f"1) {A1} = A_Sr * (1 - exp(-{lambda_Y_per_day:.4f} * T))")
print(f"2) {A2} = A_Sr * (1 - exp(-{lambda_Y_per_day:.4f} * (T + {delta_t_days})))")
print("\nDividing the two equations gives:")
print(f"{A1} / {A2} = [1 - exp(-{lambda_Y_per_day:.4f} * T)] / [1 - exp(-{lambda_Y_per_day:.4f} * T) * {k:.4f}]")

print(f"\nSolving this equation for T yields:")
print(f"The approximate time between irradiation and the first analysis is {T_days:.2f} days.")
print("\nFinal numerical answer:")
print(f"<<<{T_days:.2f}>>>")