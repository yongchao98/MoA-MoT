import math

def solve_power_optimization():
    """
    Calculates the optimal compensation power to minimize grid losses while meeting demand.
    """
    # System parameters
    P_wave_max = 5.0  # MW
    P0_wind = 10.0     # MW
    P_grid_min = 12.0  # MW
    Pc_min = -2.0      # MW
    Pc_max = 3.0       # MW
    k = 0.01         # Loss coefficient

    # The constraint on grid power is P_grid = P_bus - P_loss >= P_grid_min
    # P_bus - k * P_bus^2 >= P_grid_min
    # This can be rewritten as a quadratic inequality:
    # k * P_bus^2 - P_bus + P_grid_min <= 0
    # We solve for the roots of k*x^2 - x + P_grid_min = 0 to find the boundaries for P_bus.
    a = k
    b = -1.0
    c = P_grid_min
    
    discriminant = b**2 - 4 * a * c
    
    if discriminant < 0:
        print("Error: The system can never deliver the required minimum power.")
        return

    # The parabola opens upwards, so the valid region for P_bus is between the two roots.
    P_bus_req_lower = (-b - math.sqrt(discriminant)) / (2 * a)
    P_bus_req_upper = (-b + math.sqrt(discriminant)) / (2 * a)

    # The bounds on Pc also create a range for P_bus, since P_bus = P0_wind + P_wave + Pc
    P_bus_base = P0_wind + P_wave_max
    P_bus_from_Pc_lower = P_bus_base + Pc_min
    P_bus_from_Pc_upper = P_bus_base + Pc_max

    # The feasible region for P_bus is the intersection of the two ranges.
    feasible_P_bus_lower = max(P_bus_req_lower, P_bus_from_Pc_lower)
    feasible_P_bus_upper = min(P_bus_req_upper, P_bus_from_Pc_upper)

    if feasible_P_bus_lower > feasible_P_bus_upper:
        print("Error: There is no feasible solution that satisfies both constraints.")
        return

    # To minimize loss (k * P_bus^2), we must minimize P_bus.
    # The optimal P_bus is the lowest value in the feasible range.
    P_bus_opt = feasible_P_bus_lower
    
    # From P_bus_opt, we calculate the other optimal values.
    Pc_opt = P_bus_opt - P_bus_base
    P_loss_opt = k * P_bus_opt**2
    P_grid_opt = P_bus_opt - P_loss_opt

    print(f"The optimal compensation power required is Pc_opt = {Pc_opt:.2f} MW.")
    print("The final calculation for the total power delivered to the grid is:")
    # Print the equation with all numbers filled in as requested.
    print(f"P_grid = P_bus_opt - P_loss_opt")
    print(f"P_grid = {P_bus_opt:.2f} MW - {k} * ({P_bus_opt:.2f} MW)^2")
    print(f"P_grid = {P_bus_opt:.2f} MW - {P_loss_opt:.2f} MW")
    print(f"P_grid = {P_grid_opt:.2f} MW")

solve_power_optimization()