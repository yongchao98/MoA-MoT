import numpy as np

# --- 1. System Parameters ---
P0_wind = 10.0      # Initial wind power (MW)
P_wave_max = 5.0    # Maximum wave power (MW)
P_grid_min = 12.0   # Minimum power to grid (MW)
Pc_min = -2.0       # Minimum compensation power (MW)
Pc_max = 3.0        # Maximum compensation power (MW)
k = 0.01            # Loss coefficient
V = 1.0             # Bus voltage (p.u.)

print("### Step-by-Step Solution ###\n")

# --- 2. Formulate Power Equations ---
print("1. Power Formulations:")
print(f"Total power at the bus is P_bus = P^0_wind + P_wave.")
print(f"Power loss in the cable is P_loss = k * (P_bus / V)^2.")
print(f"Power delivered to the grid is P_grid = P_bus - P_loss.")
print("-" * 30)

# --- 3. Determine Required Bus Power ---
print("2. Determine Required Bus Power to Meet Grid Demand:")
# We need to solve P_bus - k * (P_bus/V)^2 >= P_grid_min
# This is a quadratic inequality: -k/V^2 * P_bus^2 + P_bus - P_grid_min >= 0
# Or: k/V^2 * P_bus^2 - P_bus + P_grid_min <= 0
a = k / V**2
b = -1.0
c = P_grid_min

# Find roots of the quadratic equation ax^2 + bx + c = 0
discriminant = b**2 - 4 * a * c
if discriminant < 0:
    print("Error: The system can never meet the minimum grid demand.")
else:
    P_bus_root1 = (-b - np.sqrt(discriminant)) / (2 * a)
    P_bus_root2 = (-b + np.sqrt(discriminant)) / (2 * a)

    P_bus_req_min = min(P_bus_root1, P_bus_root2)
    P_bus_req_max = max(P_bus_root1, P_bus_root2)
    
    print(f"The grid demand requires the bus power to be in the range: [{P_bus_req_min:.3f}, {P_bus_req_max:.3f}] MW.")

    # --- 4. Optimize P_wave ---
    # To minimize loss, we must minimize P_bus. So, P_bus_opt = P_bus_req_min
    P_bus_opt = P_bus_req_min
    print(f"\n3. Optimize Power Generation to Minimize Loss:")
    print(f"To minimize cable losses, we select the minimum possible bus power: P_bus_opt = {P_bus_opt:.3f} MW.")
    
    # Calculate the required wave power: P_bus = P0_wind + P_wave => P_wave = P_bus - P0_wind
    P_wave_opt = P_bus_opt - P0_wind
    
    print(f"This requires a wave power generation of P_wave_opt = {P_bus_opt:.3f} - {P0_wind:.3f} = {P_wave_opt:.3f} MW.")
    
    if P_wave_opt > P_wave_max:
        print(f"Error: The required wave power ({P_wave_opt:.3f} MW) exceeds the maximum available ({P_wave_max:.3f} MW).")
    else:
        print(f"This is feasible as it's less than the maximum wave power of {P_wave_max:.3f} MW.")
        print("-" * 30)

        # --- 5. Optimize P_c ---
        print("4. Determine Optimal Compensation Power (Pc):")
        # Pc is constrained by its given limits and the physical limit of the wave generation system
        # P_wave + Pc <= P_wave_max => Pc <= P_wave_max - P_wave
        # P_wave + Pc >= 0 => Pc >= -P_wave
        
        pc_feasible_min_1 = Pc_min
        pc_feasible_max_1 = Pc_max
        
        pc_feasible_min_2 = -P_wave_opt
        pc_feasible_max_2 = P_wave_max - P_wave_opt

        # Combine constraints
        final_pc_min = max(pc_feasible_min_1, pc_feasible_min_2)
        final_pc_max = min(pc_feasible_max_1, pc_feasible_max_2)
        
        print(f"The feasible range for Pc is [{final_pc_min:.3f}, {final_pc_max:.3f}] MW.")
        
        # We assume a secondary objective to minimize control effort, |Pc|.
        # The value in the range closest to 0 is 0 itself.
        if final_pc_min <= 0 <= final_pc_max:
            Pc_opt = 0.0
            print("Assuming a secondary objective to minimize control effort, the optimal Pc is 0.0 MW.")
        else:
            # If 0 is not in the range, choose the boundary value closest to 0.
            Pc_opt = final_pc_min if abs(final_pc_min) < abs(final_pc_max) else final_pc_max
            print(f"To minimize control effort, the optimal Pc is the boundary value {Pc_opt:.3f} MW.")

        print("-" * 30)

        # --- 6. Final Calculations and Output ---
        print("5. Final Power Delivery Calculation:")
        
        P_loss_final = k * (P_bus_opt / V)**2
        P_grid_final = P_bus_opt - P_loss_final
        
        print("\nThe final equation for power delivered to the grid is:")
        print(f"P_grid = (P^0_wind + P_wave) - k * ((P^0_wind + P_wave) / V)^2")
        print("Substituting the optimal values:")
        # We need to print each number in the equation.
        print(f"P_grid = ({P0_wind:.3f} + {P_wave_opt:.3f}) - {k:.3f} * (({P0_wind:.3f} + {P_wave_opt:.3f}) / {V:.3f})^2")
        print(f"P_grid = {P_bus_opt:.3f} - {k:.3f} * ({P_bus_opt:.3f} / {V:.3f})^2")
        print(f"P_grid = {P_bus_opt:.3f} - {P_loss_final:.3f}")
        print(f"P_grid = {P_grid_final:.3f} MW\n")

        print("### Optimal Values ###")
        print(f"Optimal Compensation Power (Pc_opt): {Pc_opt:.3f} MW")
        print(f"Total Power Delivered to Grid (P_grid): {P_grid_final:.3f} MW")

        # The final answer in the requested format
        print(f"\n<<<Optimal Compensation Power: {Pc_opt:.3f}, Total Power Delivered: {P_grid_final:.3f}>>>")
