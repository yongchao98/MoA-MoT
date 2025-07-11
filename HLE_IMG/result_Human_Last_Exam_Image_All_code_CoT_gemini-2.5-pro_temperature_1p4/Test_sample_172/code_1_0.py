import math

def solve_power_optimization():
    """
    Solves the hybrid power system optimization problem.
    """
    # Given parameters
    P_wave_max = 5  # MW
    P0_wind = 10     # MW
    P_grid_min = 12  # MW
    Pc_min = -2      # MW
    Pc_max = 3       # MW
    k = 0.01         # Loss coefficient factor
    V_bus = 1.0      # p.u.

    # We assume the wave sector operates at its maximum power for this optimization scenario.
    P_wave = P_wave_max

    # --- Step 1: Formulate the optimization problem ---
    # The total power at the system bus is P_bus = P0_wind + P_wave - Pc
    # P_bus = 10 + 5 - Pc = 15 - Pc
    # The objective is to minimize loss: P_loss = k * (P_bus / V_bus)^2 = k * P_bus^2
    # To minimize P_loss, we must minimize |P_bus|.

    # The constraint is that power delivered to the grid P_grid >= P_grid_min
    # P_grid = P_bus - P_loss = P_bus - k * P_bus^2
    # So, P_bus - k * P_bus^2 >= P_grid_min
    # k*P_bus^2 - P_bus + P_grid_min <= 0

    # --- Step 2: Solve the quadratic inequality for P_bus ---
    # We solve for the roots of a*x^2 + b*x + c = 0 where x is P_bus,
    # a = k, b = -1, c = P_grid_min
    a = k
    b = -1
    c = P_grid_min

    discriminant = b**2 - 4*a*c
    if discriminant < 0:
        print("Error: The system can never meet the minimum grid demand.")
        return

    sqrt_discriminant = math.sqrt(discriminant)
    p_bus_root1 = (-b - sqrt_discriminant) / (2*a)
    p_bus_root2 = (-b + sqrt_discriminant) / (2*a)

    # Since the parabola opens upwards (a > 0), the inequality is satisfied between the roots.
    p_bus_feasible_min_from_grid = p_bus_root1
    p_bus_feasible_max_from_grid = p_bus_root2

    # --- Step 3: Determine P_bus range from Pc limits ---
    # P_bus = 15 - Pc
    # When Pc is max (3), P_bus is min: 15 - 3 = 12
    # When Pc is min (-2), P_bus is max: 15 - (-2) = 17
    p_bus_min_from_pc = 15 - Pc_max
    p_bus_max_from_pc = 15 - Pc_min

    # --- Step 4: Find the overall feasible range for P_bus ---
    # The final feasible range is the intersection of the two ranges.
    final_feasible_min = max(p_bus_feasible_min_from_grid, p_bus_min_from_pc)
    final_feasible_max = min(p_bus_feasible_max_from_grid, p_bus_max_from_pc)

    if final_feasible_min > final_feasible_max:
        print("Error: No feasible solution found that satisfies all constraints.")
        return

    # --- Step 5: Find the optimal P_bus that minimizes loss ---
    # To minimize P_loss = k * P_bus^2, we must minimize |P_bus|.
    # Since the feasible range [final_feasible_min, final_feasible_max] is positive,
    # we choose the smallest value in the range.
    P_bus_opt = final_feasible_min

    # --- Step 6: Calculate final results ---
    # Optimal compensation power
    Pc_opt = P0_wind + P_wave - P_bus_opt
    
    # Optimal power loss
    P_loss_opt = k * (P_bus_opt**2) / (V_bus**2)
    
    # Power delivered to the grid
    P_grid_opt = P_bus_opt - P_loss_opt

    # --- Step 7: Print the results ---
    print("--- System Parameters ---")
    print(f"Initial Wind Power (P^0_wind): {P0_wind} MW")
    print(f"Wave Sector Power (P_wave): {P_wave} MW")
    print(f"Minimum Grid Delivery: {P_grid_min} MW")
    print(f"Compensation Power (P_c) Range: [{Pc_min}, {Pc_max}] MW")
    print("\n--- Optimization Results ---")
    print(f"Optimal Compensation Power (Pc_opt): {Pc_opt:.4f} MW")
    
    print("\n--- Final Power Equation Breakdown ---")
    print(f"Total power generated at the bus (P_bus):")
    print(f"P_bus = P^0_wind + P_wave - Pc_opt = {P0_wind} + {P_wave} - {Pc_opt:.4f} = {P_bus_opt:.4f} MW")
    
    print(f"\nUndersea cable power loss (P_loss):")
    print(f"P_loss = k * P_bus^2 = {k} * ({P_bus_opt:.4f})^2 = {P_loss_opt:.4f} MW")
    
    print(f"\nTotal power delivered to the onshore grid (P_grid):")
    print(f"P_grid = P_bus - P_loss = {P_bus_opt:.4f} - {P_loss_opt:.4f} = {P_grid_opt:.4f} MW")
    
    print(f"\n<<< {Pc_opt:.4f} >>>")

solve_power_optimization()