import math

# --- Given Parameters ---
P_wave_max = 5  # MW, maximum wave power
P_wind_0 = 10   # MW, initial wind power
P_grid_min = 12 # MW, minimum power to be delivered to the grid
Pc_min = -2     # MW, minimum compensation power
Pc_max = 3      # MW, maximum compensation power
k = 0.01        # Proportionality constant for power loss (in 1/MW)
V_bus = 1       # p.u., bus voltage

# --- Step 1: Assume P_wave is at its maximum for this optimization scenario ---
P_wave = P_wave_max

print("Problem Analysis & Formulation:")
print("---------------------------------")
print(f"Initial Wind Power (P_wind^0): {P_wind_0} MW")
print(f"Wave Power (P_wave): {P_wave} MW (Assumed at maximum)")
print(f"Grid Demand (P_grid_min): {P_grid_min} MW")
print(f"Compensation Power (Pc) Range: [{Pc_min}, {Pc_max}] MW")
print(f"Loss coefficient (k): {k} 1/MW")
print("\nTotal power at the bus is formulated as: P_bus = P_wind^0 + P_wave + Pc")
print("Power loss is: P_loss = k * P_bus^2")
print("Power to grid is: P_grid = P_bus - P_loss")
print("---------------------------------\n")

# --- Step 2: Solve the grid demand constraint to find the required P_bus range ---
# The constraint is P_bus - k*P_bus^2 >= P_grid_min
# This rearranges to a quadratic inequality: k*P_bus^2 - P_bus + P_grid_min <= 0
a = k
b = -1
c = P_grid_min

# Calculate the discriminant
delta = b**2 - 4 * a * c

if delta < 0:
    print("The grid demand cannot be met at any power level.")
else:
    # Find the roots of the quadratic equation k*x^2 - x + c = 0
    # These roots define the boundaries where P_grid = P_grid_min
    sqrt_delta = math.sqrt(delta)
    P_bus_req1 = (-b - sqrt_delta) / (2 * a)
    P_bus_req2 = (-b + sqrt_delta) / (2 * a)
    
    # To satisfy the inequality, P_bus must be between the two roots.
    P_bus_req_min = min(P_bus_req1, P_bus_req2)
    
    # --- Step 3: Determine the optimal bus power ---
    # To minimize loss (k*P_bus^2), we must choose the minimum possible P_bus.
    P_bus_opt = P_bus_req_min

    # --- Step 4: Calculate the optimal compensation power ---
    # P_bus_opt = P_wind_0 + P_wave + Pc_opt
    Pc_opt = P_bus_opt - P_wind_0 - P_wave

    print("Calculation Results:")
    print("---------------------------------")
    
    # Check if the optimal Pc is within its allowed physical range
    if Pc_min <= Pc_opt <= Pc_max:
        print(f"To minimize cable losses while meeting grid demand, the bus power must be: {P_bus_opt:.2f} MW.")
        print(f"The optimal compensation power required is: Pc_opt = {Pc_opt:.2f} MW.")
        
        # --- Step 5: Compute the final power delivered and loss ---
        P_loss_opt = k * (P_bus_opt**2)
        P_grid_opt = P_bus_opt - P_loss_opt

        print("\n--- Final Power Delivery Calculation ---")
        # Final equation with numbers
        print(f"Total Power Delivered (P_grid) = P_bus_opt - k * (P_bus_opt)^2")
        print(f"P_grid = {P_bus_opt:.2f} MW - {k:.2f} * ({P_bus_opt:.2f})^2 MW")
        print(f"P_grid = {P_bus_opt:.2f} MW - {P_loss_opt:.2f} MW = {P_grid_opt:.2f} MW")

    else:
        print(f"The calculated optimal Pc ({Pc_opt:.2f} MW) is outside the allowed range [{Pc_min}, {Pc_max}] MW.")
        print("This means the system must operate at a non-optimal point to respect component limits.")
        # Logic for this case would follow, but based on problem values, it's not needed.
