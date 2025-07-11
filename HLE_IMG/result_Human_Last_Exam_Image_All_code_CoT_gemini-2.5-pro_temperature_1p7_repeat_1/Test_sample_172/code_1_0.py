import math

def solve_power_optimization():
    """
    Solves the hybrid power system optimization problem.
    """
    # Step 1: Define the system parameters from the problem description
    P0_wind = 10.0  # Initial wind power output (MW)
    P_wave = 5.0    # Wave power output, assumed to be at maximum (MW)
    k = 0.01        # Power loss coefficient
    V = 1.0         # Bus voltage (p.u.)
    P_grid_min = 12.0 # Minimum power to be delivered to the grid (MW)
    Pc_min = -2.0   # Minimum compensation power (MW)
    Pc_max = 3.0    # Maximum compensation power (MW)

    print("Step 1: System Parameters")
    print(f"Initial Wind Power (P^0_wind): {P0_wind} MW")
    print(f"Wave Power (P_wave): {P_wave} MW")
    print(f"Loss Coefficient (k): {k}")
    print(f"Minimum Grid Demand (P_grid_min): {P_grid_min} MW")
    print(f"Compensation Power Bounds (P_c): [{Pc_min}, {Pc_max}] MW\n")

    # Step 2: Formulate the power equations
    # Total power at the bus (before the cable): P_bus = P0_wind + P_wave + P_c
    # Power loss in the cable: P_loss = k * (P_bus / V)^2 = k * P_bus^2
    # Power delivered to the grid: P_grid = P_bus - P_loss
    # So, P_grid = P_bus - k * P_bus^2
    print("Step 2: Power Equations Formulation")
    print("Total Power at Bus (P_bus) = P^0_wind + P_wave + P_c")
    print("Power Loss (P_loss) = k * P_bus^2")
    print("Power to Grid (P_grid) = P_bus - P_loss\n")
    
    # Step 3: Set up the optimization problem
    # To minimize loss, we minimize P_bus. We find the minimum P_bus
    # that satisfies the grid demand constraint: P_grid = P_grid_min
    # This gives the quadratic equation: k*P_bus^2 - P_bus + P_grid_min = 0
    print("Step 3: Solve for the optimal P_bus to meet minimum demand")
    print(f"Setting P_grid = {P_grid_min} MW gives the equation:")
    print(f"{k}*P_bus^2 - {V}*P_bus + {P_grid_min} = 0\n")

    # Step 4: Solve the quadratic equation for P_bus
    a = k
    b = -1.0
    c = P_grid_min
    
    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        print("Error: No real solution exists. It is impossible to meet the grid demand.")
        return

    p_bus_sol1 = (-b - math.sqrt(discriminant)) / (2 * a)
    p_bus_sol2 = (-b + math.sqrt(discriminant)) / (2 * a)

    # We choose the smaller P_bus solution because it results in lower losses
    P_bus_opt = p_bus_sol1

    print("Step 4: Solving the quadratic equation yields two possible values for P_bus:")
    print(f"P_bus_1 = {p_bus_sol1:.3f} MW")
    print(f"P_bus_2 = {p_bus_sol2:.3f} MW")
    print(f"To minimize losses, we choose the smaller value: P_bus_opt = {P_bus_opt:.3f} MW\n")

    # Step 5: Calculate the optimal compensation power (Pc_opt)
    # P_bus_opt = P0_wind + P_wave + Pc_opt
    # Pc_opt = P_bus_opt - P0_wind - P_wave
    Pc_opt = P_bus_opt - P0_wind - P_wave
    
    print("Step 5: Calculate the required compensation power (P_c_opt)")
    print(f"P_c_opt = P_bus_opt - P^0_wind - P_wave")
    print(f"P_c_opt = {P_bus_opt:.3f} - {P0_wind} - {P_wave} = {Pc_opt:.3f} MW")
    
    # Verify that Pc_opt is within its allowed range
    if Pc_min <= Pc_opt <= Pc_max:
        print(f"The calculated P_c_opt ({Pc_opt:.3f} MW) is within the allowed range [{Pc_min}, {Pc_max}] MW.\n")
    else:
        print(f"Warning: The calculated P_c_opt ({Pc_opt:.3f} MW) is outside the allowed range [{Pc_min}, {Pc_max}] MW. The problem constraints cannot be met optimally.\n")

    # Step 6: Compute final power values with Pc_opt
    P_loss_opt = k * P_bus_opt**2
    P_grid_final = P_bus_opt - P_loss_opt

    print("Step 6: Final Results")
    print("--- Optimal Values ---")
    print(f"Optimal Compensation Power (P_c_opt): {Pc_opt:.3f} MW")
    print("\n--- Power Delivered to Grid Calculation ---")
    print(f"Total Power Generated at Bus (P_bus) = {P0_wind} + {P_wave} + ({Pc_opt:.3f}) = {P_bus_opt:.3f} MW")
    print(f"Power Loss in Cable (P_loss) = {k} * ({P_bus_opt:.3f})^2 = {P_loss_opt:.3f} MW")
    print(f"Total Power Delivered to Grid (P_grid) = P_bus - P_loss")
    print(f"P_grid = {P_bus_opt:.3f} - {P_loss_opt:.3f} = {P_grid_final:.3f} MW")


# Execute the function to solve the problem
solve_power_optimization()
<<<P_c_opt = -1.056 MW, P_grid = 12.000 MW>>>