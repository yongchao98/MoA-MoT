import math

# Step 1: Define System Variables and Constants
P_wave_max = 5  # MW
P0_wind = 10    # MW
P_grid_min = 12 # MW
Pc_min = -2     # MW
Pc_max = 3      # MW
k = 0.01        # Proportionality constant for power loss

print("--- System Parameters ---")
print(f"Initial Wind Power (P0_wind): {P0_wind} MW")
print(f"Maximum Wave Power (P_wave_max): {P_wave_max} MW")
print(f"Minimum Grid Demand (P_grid_min): {P_grid_min} MW")
print(f"Compensation Power Range (Pc): [{Pc_min}, {Pc_max}] MW")
print(f"Loss Coefficient (k): {k}")
print("-" * 27 + "\n")

# Step 2: Formulate Power Equations (as explained in the plan)
# P_bus = P0_wind + P_wave + Pc
# P_loss = k * P_bus^2
# P_grid = P_bus - P_loss

# Step 3: Solve the Constraint for P_bus
# We need P_grid >= P_grid_min, which means:
# P_bus - k * P_bus^2 >= P_grid_min
# This can be rewritten as a standard quadratic inequality:
# k * P_bus^2 - P_bus + P_grid_min <= 0

# To find the roots of k*x^2 - x + P_grid_min = 0 where x = P_bus
a = k
b = -1
c = P_grid_min

# Using the quadratic formula: x = [-b +/- sqrt(b^2 - 4ac)] / 2a
discriminant = b**2 - 4 * a * c

if discriminant < 0:
    print("The grid demand cannot be met.")
else:
    P_bus_root1 = (-b - math.sqrt(discriminant)) / (2 * a)
    P_bus_root2 = (-b + math.sqrt(discriminant)) / (2 * a)

    # The quadratic opens upwards, so the inequality is satisfied between the roots.
    P_bus_req_min = P_bus_root1
    P_bus_req_max = P_bus_root2
    print(f"--- Power Constraint Analysis ---")
    print(f"To meet grid demand, P_bus must be in the range: [{P_bus_req_min:.3f}, {P_bus_req_max:.3f}] MW")

    # Step 4: Incorporate Pc Constraints
    # Calculate the range of P_bus achievable with the given Pc limits
    P_base = P0_wind + P_wave_max
    P_bus_achievable_min = P_base + Pc_min
    P_bus_achievable_max = P_base + Pc_max
    print(f"The achievable P_bus range with Pc is: [{P_bus_achievable_min:.3f}, {P_bus_achievable_max:.3f}] MW")

    # Step 5: Find the Optimal P_bus
    # The overall valid range for P_bus is the intersection of the required and achievable ranges.
    valid_P_bus_min = max(P_bus_req_min, P_bus_achievable_min)
    valid_P_bus_max = min(P_bus_req_max, P_bus_achievable_max)

    # To minimize loss (P_loss = k * P_bus^2), we must minimize P_bus.
    # So, the optimal P_bus is the minimum value in the valid range.
    P_bus_opt = valid_P_bus_min
    
    print(f"The combined valid range for P_bus is: [{valid_P_bus_min:.3f}, {valid_P_bus_max:.3f}] MW")
    print("-" * 31 + "\n")

    # Step 6: Calculate Final Results
    # Calculate the optimal compensation power
    Pc_opt = P_bus_opt - P_base
    
    # Calculate the power loss and grid power at the optimal point
    P_loss_opt = k * P_bus_opt**2
    P_grid_opt = P_bus_opt - P_loss_opt

    print("--- Optimal Power Configuration ---")
    print(f"Optimal Compensation Power (Pc_opt): {Pc_opt:.3f} MW")
    print(f"This results in an optimal bus power (P_bus_opt) of: {P_bus_opt:.3f} MW")
    
    print("\n--- Final Power Delivery ---")
    print(f"Power Loss at Optimum (P_loss_opt): {P_loss_opt:.3f} MW")
    print(f"Total Power Delivered to Grid (P_grid_opt): {P_grid_opt:.3f} MW")

    print("\nFinal Equation (P_grid = P_bus - P_loss):")
    print(f"{P_grid_opt:.3f} MW = {P_bus_opt:.3f} MW - {P_loss_opt:.3f} MW")

    # The final numerical answer for Pc_opt
    print(f"\n<<<{Pc_opt:.3f}>>>")