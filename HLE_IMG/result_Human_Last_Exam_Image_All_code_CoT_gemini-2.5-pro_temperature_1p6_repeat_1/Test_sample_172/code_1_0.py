import math

def solve_power_system_optimization():
    """
    Solves for the optimal compensation power to meet grid demand
    while minimizing undersea cable losses.
    """
    # 1. Define constants from the problem statement
    P_wave_max = 5  # MW
    P0_wind = 10  # MW
    P_grid_demand = 12  # MW
    Pc_min = -2  # MW
    Pc_max = 3  # MW
    k = 0.01  # Loss coefficient
    V_bus = 1  # p.u.

    # Assume the system operates at the maximum wave power given
    P_wave = P_wave_max
    
    # Based on the analysis, the total power on the bus is modeled as:
    # P_bus = P0_wind + P_wave - P_c
    P_base_gen = P0_wind + P_wave
    
    print("Step 1: Formulate the power equations.")
    print(f"Total available power from sources P_avail = P0_wind + P_wave = {P0_wind} + {P_wave} = {P_base_gen} MW.")
    print("Based on the problem description, we model the compensated bus power as P_bus = P_avail - P_c.")
    print("The power delivered to the grid is P_grid = P_bus - P_loss, where P_loss = k * P_bus^2.")
    print(f"Thus, P_grid = (P_avail - P_c) - k * (P_avail - P_c)^2 >= {P_grid_demand} MW.\n")

    print("Step 2: Solve the grid demand constraint to find the required P_bus.")
    # We need to solve the quadratic equation P_bus - k * P_bus^2 = P_grid_demand
    # k*P_bus^2 - P_bus + P_grid_demand = 0
    a = k
    b = -1
    c = P_grid_demand
    
    # Using the quadratic formula: P_bus = [-b ± sqrt(b^2 - 4ac)] / 2a
    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        print("Error: The system cannot meet the grid demand under any condition.")
        return

    sqrt_discriminant = math.sqrt(discriminant)
    
    P_bus_sol1 = (-b - sqrt_discriminant) / (2 * a)
    P_bus_sol2 = (-b + sqrt_discriminant) / (2 * a)

    # For the inequality P_bus - k*P_bus^2 >= P_grid_demand, P_bus must be between the two solutions.
    P_bus_min_req = P_bus_sol1
    P_bus_max_req = P_bus_sol2
    print(f"To meet grid demand, P_bus must be in the range: [{P_bus_min_req:.3f}, {P_bus_max_req:.3f}] MW.\n")

    print("Step 3: Determine the optimal compensation power P_c.")
    # The objective is to minimize P_loss = k * P_bus^2. This means minimizing P_bus.
    # To minimize P_bus = P_base_gen - P_c, we must maximize P_c.
    # From the constraint P_bus >= P_bus_min_req, we have:
    # P_base_gen - P_c >= P_bus_min_req
    # P_c <= P_base_gen - P_bus_min_req
    Pc_upper_bound_from_grid = P_base_gen - P_bus_min_req
    
    # The final valid range for P_c is the intersection of its original bounds and the new derived bound.
    # [Pc_min, Pc_max] intersected with (-inf, Pc_upper_bound_from_grid]
    final_Pc_min = Pc_min
    final_Pc_max = min(Pc_max, Pc_upper_bound_from_grid)
    
    print(f"The grid demand requires P_c <= {Pc_upper_bound_from_grid:.3f} MW.")
    print(f"Combined with the initial bounds, the valid range for P_c is [{final_Pc_min}, {final_Pc_max:.3f}] MW.")
    
    # To maximize P_c, we choose the highest value in the valid range.
    Pc_opt = final_Pc_max
    print(f"To minimize loss, we must maximize P_c. Therefore, the optimal compensation power is P_c_opt = {Pc_opt:.3f} MW.\n")
    
    print("Step 4: Calculate the final power values with optimal P_c.")
    # Calculate final system state with optimal P_c
    P_bus_opt = P_base_gen - Pc_opt
    P_loss_opt = k * P_bus_opt**2
    P_grid_final = P_bus_opt - P_loss_opt
    
    print(f"Optimal Bus Power (P_bus) = {P_base_gen} - {Pc_opt:.3f} = {P_bus_opt:.3f} MW")
    print(f"Minimized Cable Loss (P_loss) = {k} * ({P_bus_opt:.3f})^2 = {P_loss_opt:.3f} MW")
    print("\n--- Final Results ---")
    print(f"The optimal compensation power required is P_c_opt = {Pc_opt:.3f} MW.")
    print("The final equation for power delivered to the grid is:")
    # Using each number in the final equation format as requested
    print(f"{P_bus_opt:.3f} (P_bus) - {P_loss_opt:.3f} (P_loss) = {P_grid_final:.3f} (P_grid) MW")
    

solve_power_system_optimization()
<<<1.056>>>