import math

def solve_power_optimization():
    """
    Calculates the optimal compensation power and grid power for the hybrid offshore system.
    """
    # System parameters
    P_wave_max = 5.0  # MW
    P0_wind = 10.0   # MW
    P_grid_min = 12.0  # MW
    Pc_min = -2.0      # MW
    Pc_max = 3.0       # MW
    k = 0.01

    # --- Step 1: Find the feasible range for P_bus ---
    # We solve the quadratic inequality: k*P_bus^2 - P_bus + P_grid_min <= 0
    # The corresponding equation is k*x^2 - x + P_grid_min = 0
    a = k
    b = -1
    c = P_grid_min
    
    # Calculate the discriminant
    discriminant = b**2 - 4*a*c

    if discriminant < 0:
        print("Error: The system cannot meet the minimum grid demand under any condition.")
        return

    # Find the roots of the quadratic equation
    sqrt_discriminant = math.sqrt(discriminant)
    p_bus_root1 = (-b - sqrt_discriminant) / (2*a)
    p_bus_root2 = (-b + sqrt_discriminant) / (2*a)
    
    # The valid range for P_bus is between the two roots
    p_bus_min_req = min(p_bus_root1, p_bus_root2)
    
    # --- Step 2: Determine optimal P_bus and P_wave ---
    # To minimize loss (k*P_bus^2), we choose the minimum valid P_bus
    P_bus_opt = p_bus_min_req
    
    # Required wave power is P_wave_opt = P_bus_opt - P0_wind
    P_wave_opt = P_bus_opt - P0_wind

    if not (0 <= P_wave_opt <= P_wave_max):
        print(f"Error: Optimal wave power ({P_wave_opt:.3f} MW) is outside the feasible range [0, {P_wave_max} MW].")
        return

    # --- Step 3: Determine optimal Pc ---
    # Assume a secondary objective: run the wave converters at max rated power output
    # P_wave_out = P_wave_opt + Pc_opt = P_wave_max
    Pc_opt = P_wave_max - P_wave_opt

    if not (Pc_min <= Pc_opt <= Pc_max):
         print(f"Error: Calculated Pc_opt ({Pc_opt:.3f} MW) is outside the feasible range [{Pc_min}, {Pc_max} MW].")
         return

    # --- Step 4: Calculate final results and print ---
    P_loss_opt = k * P_bus_opt**2
    P_grid_opt = P_bus_opt - P_loss_opt

    print(f"The optimal compensation power required is: Pc_opt = {Pc_opt:.3f} MW\n")
    print("--- Calculation of Total Power Delivered to Grid ---")
    print(f"P_grid = (P^0_wind + P_wave_opt) - k * (P^0_wind + P_wave_opt)^2")
    print(f"P_grid = ({P0_wind:.1f} + {P_wave_opt:.3f}) - {k} * ({P0_wind:.1f} + {P_wave_opt:.3f})^2")
    print(f"P_grid = {P_bus_opt:.3f} - {k} * ({P_bus_opt:.3f})^2")
    print(f"P_grid = {P_bus_opt:.3f} - {P_loss_opt:.3f}")
    print(f"P_grid = {P_grid_opt:.3f} MW")

    # Final answer in the required format
    print(f"\n<<<{Pc_opt:.3f}>>>")

# Execute the function
solve_power_optimization()