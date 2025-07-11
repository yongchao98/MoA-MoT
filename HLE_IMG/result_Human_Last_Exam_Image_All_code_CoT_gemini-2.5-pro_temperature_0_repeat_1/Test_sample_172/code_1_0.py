import math

def solve_power_optimization():
    """
    Solves for the optimal compensation power to minimize cable losses
    while meeting grid demand.
    """
    # Given parameters
    P_wave_max = 5.0  # MW
    P0_wind = 10.0    # MW
    P_grid_min = 12.0 # MW
    Pc_min = -2.0     # MW
    Pc_max = 3.0      # MW
    k = 0.01          # Loss coefficient
    V = 1.0           # p.u. voltage

    # The total power at the bus is P_bus = P0_wind + P_wave + Pc
    # P_bus = 10 + 5 + Pc = 15 + Pc
    # The power delivered to the grid is P_grid = P_bus - k * P_bus^2
    # We need P_grid >= 12, which means:
    # P_bus - k * P_bus^2 >= 12
    # k * P_bus^2 - P_bus + 12 <= 0
    # 0.01 * P_bus^2 - P_bus + 12 <= 0

    # Solve the quadratic equation 0.01*x^2 - x + 12 = 0 for x = P_bus
    a = k
    b = -1
    c = P_grid_min
    
    discriminant = b**2 - 4*a*c
    if discriminant < 0:
        print("The grid demand can never be met.")
        return

    sqrt_discriminant = math.sqrt(discriminant)
    P_bus_root1 = (-b - sqrt_discriminant) / (2*a)
    P_bus_root2 = (-b + sqrt_discriminant) / (2*a)

    # The quadratic opens upwards, so the solution is between the roots.
    P_bus_feasible_min = P_bus_root1
    P_bus_feasible_max = P_bus_root2

    # Convert the feasible P_bus range to a feasible Pc range
    # Pc = P_bus - (P0_wind + P_wave_max)
    Pc_feasible_min_from_grid = P_bus_feasible_min - (P0_wind + P_wave_max)
    Pc_feasible_max_from_grid = P_bus_feasible_max - (P0_wind + P_wave_max)

    # Find the intersection of the feasible Pc range and the given Pc limits
    final_Pc_min = max(Pc_min, Pc_feasible_min_from_grid)
    final_Pc_max = min(Pc_max, Pc_feasible_max_from_grid)

    if final_Pc_min > final_Pc_max:
        print("No feasible solution for Pc exists that satisfies all constraints.")
        return

    # To minimize loss P_loss = k * P_bus^2, we need to minimize P_bus.
    # P_bus = 15 + Pc. To minimize P_bus, we must choose the minimum possible Pc.
    Pc_opt = final_Pc_min

    # Calculate the final values using the optimal Pc
    P_bus_opt = P0_wind + P_wave_max + Pc_opt
    P_loss_opt = k * P_bus_opt**2
    P_grid_opt = P_bus_opt - P_loss_opt

    # Print the formulation of total power delivered to the grid
    print("Formulation for total power delivered to the grid (P_grid):")
    print("P_grid = (P^0_wind + P_wave + P_c) - k * (P^0_wind + P_wave + P_c)^2")
    print(f"P_grid = ({P0_wind} + {P_wave_max} + P_c) - {k} * ({P0_wind} + {P_wave_max} + P_c)^2\n")

    # Print the results
    print("--- Optimal Power Configuration ---")
    print(f"Optimal Compensation Power (Pc_opt): {Pc_opt:.4f} MW")
    print(f"Total Power at System Bus (P_bus): {P_bus_opt:.4f} MW")
    print(f"Undersea Cable Power Loss (P_loss): {P_loss_opt:.4f} MW")
    print(f"Total Power Delivered to Grid (P_grid): {P_grid_opt:.4f} MW\n")
    
    # Print the final equation with numbers
    print("Final power delivery equation:")
    print(f"P_grid = P_bus - P_loss")
    print(f"{P_grid_opt:.4f} MW = {P_bus_opt:.4f} MW - {P_loss_opt:.4f} MW")


solve_power_optimization()
<<< -1.0555 >>>