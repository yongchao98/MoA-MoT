import math

# --- Given Parameters ---
P_wind_0 = 10  # MW
P_wave_max = 5  # MW
P_demand_min = 12  # MW
Pc_min = -2  # MW
Pc_max = 3  # MW
k = 0.01  # Loss coefficient
V = 1  # p.u.

# For this problem, we assume the wave sector operates at its maximum power.
P_wave = P_wave_max

# --- 1. Formulate the general equations ---
print("Step 1: Formulating the Power Equations")
print("Total power at the bus (P_bus) = P_wind_0 + P_wave + P_c")
print(f"Total power delivered to the grid (P_grid) = P_bus - k * (P_bus / V)^2\n")

# --- 2. Solve the Grid Demand Constraint ---
# The constraint is P_bus - k * P_bus^2 >= P_demand_min
# Rearranging gives a quadratic inequality: k*P_bus^2 - P_bus + P_demand_min <= 0
# We solve for the roots of a*x^2 + b*x + c = 0 where x = P_bus
a = k
b = -1
c = P_demand_min

# Using the quadratic formula: x = [-b +/- sqrt(b^2 - 4ac)] / 2a
discriminant = b**2 - 4 * a * c
if discriminant < 0:
    print("Error: The grid demand cannot be met.")
else:
    sqrt_discriminant = math.sqrt(discriminant)
    P_bus_demand_min = (-b - sqrt_discriminant) / (2 * a)
    P_bus_demand_max = (-b + sqrt_discriminant) / (2 * a)

    print("Step 2: Solving the Grid Demand Constraint")
    print(f"The inequality {k}*P_bus^2 - {1}*P_bus + {P_demand_min} <= 0 must be satisfied.")
    print(f"This requires P_bus to be between {P_bus_demand_min:.3f} MW and {P_bus_demand_max:.3f} MW.\n")

    # --- 3. Determine the Supply Constraint ---
    P_bus_base = P_wind_0 + P_wave
    P_bus_supply_min = P_bus_base + Pc_min
    P_bus_supply_max = P_bus_base + Pc_max

    print("Step 3: Determining the Supply Constraint")
    print(f"With P_c ranging from {Pc_min} to {Pc_max} MW, the possible P_bus range is [{P_bus_supply_min:.3f}, {P_bus_supply_max:.3f}] MW.\n")

    # --- 4. Find the Optimal Bus Power (P_bus_opt) ---
    # The feasible range for P_bus is the intersection of the demand and supply ranges.
    # We want to minimize P_bus, so we take the maximum of the lower bounds.
    P_bus_opt = max(P_bus_demand_min, P_bus_supply_min)
    
    # Check if a solution exists
    if P_bus_opt > min(P_bus_demand_max, P_bus_supply_max):
        print("Error: No feasible solution exists that satisfies both constraints.")
    else:
        print("Step 4: Finding the Optimal Solution")
        print(f"To minimize losses, we must minimize P_bus.")
        print(f"The optimal bus power P_bus_opt is the minimum value that satisfies both constraints: {P_bus_opt:.3f} MW.\n")

        # --- 5. Calculate Final Values ---
        # Optimal compensation power
        Pc_opt = P_bus_opt - P_bus_base

        # Total power delivered to the grid
        P_loss_opt = k * (P_bus_opt / V)**2
        P_grid_opt = P_bus_opt - P_loss_opt
        
        print("--- Final Results ---")
        print(f"Optimal Compensation Power (Pc_opt): {Pc_opt:.3f} MW")
        print("\nFinal Equation for Optimal Bus Power:")
        print(f"P_bus_opt = P_wind_0 + P_wave + Pc_opt")
        print(f"{P_bus_opt:.3f} MW = {P_wind_0} MW + {P_wave} MW + ({Pc_opt:.3f} MW)")
        
        print("\nFinal Equation for Power Delivered to Grid:")
        print(f"P_grid_opt = P_bus_opt - k * P_bus_opt^2")
        print(f"{P_grid_opt:.3f} MW = {P_bus_opt:.3f} MW - {k} * ({P_bus_opt:.3f} MW)^2 = {P_bus_opt:.3f} MW - {P_loss_opt:.3f} MW")

        # --- Answer for automated checker ---
        print(f"\n<<<{Pc_opt:.3f}>>>")