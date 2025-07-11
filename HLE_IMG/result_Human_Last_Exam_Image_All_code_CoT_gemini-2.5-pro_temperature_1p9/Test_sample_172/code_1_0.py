import math

def solve_power_optimization():
    """
    Solves the hybrid power system optimization problem to determine the optimal
    compensation power and the final power delivered to the grid, while providing
    a step-by-step calculation.
    """
    # Step 0: Define constants from the problem description
    P_wave_max = 5.0  # MW
    P_wind_0 = 10.0   # MW
    P_grid_min = 12.0 # MW
    Pc_min = -2.0     # MW
    Pc_max = 3.0      # MW
    k = 0.01          # 1/MW (loss coefficient)

    # --- Calculations ---
    
    # 1. Determine the required P_bus from the grid demand constraint.
    # The constraint is P_grid >= 12, which means: P_bus - k*P_bus^2 >= 12
    # This rearranges to a standard quadratic inequality: k*P_bus^2 - P_bus + 12 <= 0
    a, b, c = k, -1.0, P_grid_min
    # Find the roots of ax^2 + bx + c = 0 to determine the valid range for P_bus
    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        print("Error: The system cannot meet the minimum grid demand under any condition.")
        return

    sqrt_discriminant = math.sqrt(discriminant)
    # The lower root of the quadratic equation gives the minimum required bus power
    p_bus_demand_min = (-b - sqrt_discriminant) / (2 * a)
    
    # 2. Determine the possible P_bus from the compensation power (Pc) limits.
    base_power = P_wind_0 + P_wave_max
    # The minimum P_bus achievable via control is when Pc is at its minimum
    p_bus_control_min = base_power + Pc_min
    
    # 3. The optimal P_bus is the minimum value that satisfies both constraints.
    P_bus_opt = max(p_bus_demand_min, p_bus_control_min)
    
    # 4. Calculate all final values based on the optimal P_bus.
    P_c_opt = P_bus_opt - base_power
    P_loss_opt = k * P_bus_opt**2
    P_grid_final = P_bus_opt - P_loss_opt

    # --- Step-by-step Output ---
    
    print("Step 1: Determine the Optimal Bus Power (P_bus_opt)\n")
    print(f"The system must deliver at least {P_grid_min} MW to the grid.")
    print(f"This is expressed by the equation: P_grid = P_bus - {k} * P_bus^2 >= {P_grid_min}")
    print(f"Solving this inequality shows that the bus power (P_bus) must be at least {p_bus_demand_min:.3f} MW.\n")

    print(f"The compensation power (Pc) is limited to [{Pc_min}, {Pc_max}] MW.")
    print(f"This limits the total bus power to a range based on P_bus = ({P_wind_0} + {P_wave_max}) + Pc.")
    print(f"The minimum possible bus power based on control is {base_power} + {Pc_min} = {p_bus_control_min:.1f} MW.\n")
    
    print("To minimize losses, we need the smallest P_bus that satisfies BOTH conditions.")
    print(f"Therefore, P_bus_opt = max({p_bus_demand_min:.3f}, {p_bus_control_min:.1f}) = {P_bus_opt:.3f} MW.\n")
    print("-" * 60)

    print("Step 2: Calculate the Optimal Compensation Power (Pc_opt)\n")
    print("The optimal compensation power is found using the P_bus_opt value:")
    print("   Pc_opt = P_bus_opt - (P_wind_0 + P_wave)")
    print(f"   Pc_opt = {P_bus_opt:.3f} MW - ({P_wind_0} MW + {P_wave_max} MW) = {P_c_opt:.3f} MW\n")
    print("-" * 60)

    print("Step 3: Calculate the Final Power Delivered to the Grid (P_grid_final)\n")
    print("First, calculate the power loss with the optimal bus power:")
    print("   P_loss_opt = k * (P_bus_opt)^2")
    print(f"   P_loss_opt = {k} * ({P_bus_opt:.3f})^2 = {P_loss_opt:.3f} MW\n")
    
    print("Finally, calculate the total power delivered to the grid:")
    print("   P_grid_final = P_bus_opt - P_loss_opt")
    print(f"   P_grid_final = {P_bus_opt:.3f} MW - {P_loss_opt:.3f} MW = {P_grid_final:.3f} MW\n")

    global final_answer
    final_answer = P_c_opt


# Run the solver
solve_power_optimization()
# Return the optimal compensation power as the final answer
print(f"<<<{final_answer:.3f}>>>")