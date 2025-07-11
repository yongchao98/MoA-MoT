import math

def solve_power_optimization():
    """
    Calculates the optimal compensation power (Pc) and the total power
    delivered to the grid for a hybrid offshore power generation system.
    """
    # --- Given Parameters ---
    P_wave_max = 5.0  # MW
    P_wind_0 = 10.0   # MW
    P_grid_req = 12.0 # MW, minimum power to be delivered to the grid
    Pc_min = -2.0     # MW
    Pc_max = 3.0      # MW
    k = 0.01          # Loss coefficient (in MW^-1, assuming V=1 p.u. and P is in MW)

    # Assume wave sector operates at max power to contribute to grid demand
    P_wave = P_wave_max

    # --- Step 1: Define relationship between P_bus and P_grid ---
    # The objective is to minimize loss (P_loss = k * P_bus^2), which means minimizing P_bus.
    # The main constraint is that power delivered to the grid must meet the requirement:
    # P_grid = P_bus - P_loss >= P_grid_req
    # P_bus - k * P_bus^2 >= P_grid_req
    # Rearranging gives the quadratic inequality: k*P_bus^2 - P_bus + P_grid_req <= 0

    # --- Step 2: Solve the quadratic equation to find the feasible range for P_bus ---
    # We solve a*x^2 + b*x + c = 0 where x = P_bus
    a = k
    b = -1
    c = P_grid_req

    # Calculate the discriminant
    discriminant = b**2 - 4 * a * c

    if discriminant < 0:
        print("Error: The system cannot meet the minimum grid demand under any condition.")
        return

    # The inequality holds between the two roots of the quadratic equation.
    # The lower root represents the minimum P_bus required to deliver P_grid_req.
    P_bus_min_required = (-b - math.sqrt(discriminant)) / (2 * a)

    # --- Step 3: Determine the required Pc and find the optimal Pc ---
    # Calculate the Pc needed to achieve this minimum required P_bus.
    # P_bus_min_required = P_wind_0 + P_wave + Pc_required
    Pc_required = P_bus_min_required - P_wind_0 - P_wave

    # The optimal Pc is the minimum possible value that satisfies all constraints.
    # The feasible range for Pc is [max(Pc_min, Pc_required), Pc_max].
    # To minimize P_bus, we must select the lowest value in this range.
    Pc_opt = max(Pc_min, Pc_required)
    
    if Pc_opt > Pc_max:
        print(f"Error: The optimal required Pc ({Pc_opt:.3f} MW) is outside the controller's maximum limit of {Pc_max} MW.")
        print("The system cannot satisfy the grid demand.")
        return

    # --- Step 4: Calculate final power values using Pc_opt ---
    P_bus_opt = P_wind_0 + P_wave + Pc_opt
    P_loss_opt = k * (P_bus_opt**2)
    P_grid_final = P_bus_opt - P_loss_opt

    # --- Final Output ---
    print("--- Power System Optimization Calculation ---\n")
    print("1. Formulation of Total Power Delivered to the Grid:")
    print("The total power delivered (P_grid) is the bus power (P_bus) minus cable losses (P_loss).")
    print("P_grid = P_bus - P_loss")
    print("P_grid = (P^0_wind + P_wave + P_c) - k * (P^0_wind + P_wave + P_c)^2\n")

    print(f"2. Optimal Compensation Power (Pc_opt):")
    print(f"To satisfy the grid demand of {P_grid_req} MW while minimizing transmission losses,")
    print(f"the minimum required power at the system bus is {P_bus_min_required:.3f} MW.")
    print(f"The compensation power needed to achieve this is:")
    print(f"Pc_required = {P_bus_min_required:.3f} MW - {P_wind_0} MW - {P_wave} MW = {Pc_required:.3f} MW.")
    print(f"This value is within the allowed range [{Pc_min}, {Pc_max}] MW, so it is the optimal value.")
    print(f"Optimal Compensation Power, Pc_opt = {Pc_opt:.3f} MW\n")

    print("3. Total Power Delivered to the Grid (with optimal compensation):")
    print("Using Pc_opt, the final equation for the power delivered to the grid is:")
    # Print the equation with all the numbers substituted
    print(f"P_grid = ({P_wind_0} + {P_wave} + ({Pc_opt:.3f})) - {k} * ({P_wind_0} + {P_wave} + ({Pc_opt:.3f}))^2")
    print(f"P_grid = {P_bus_opt:.3f} - {k} * ({P_bus_opt:.3f})^2")
    print(f"P_grid = {P_bus_opt:.3f} - {k} * {P_bus_opt**2:.3f}")
    print(f"P_grid = {P_bus_opt:.3f} - {P_loss_opt:.3f}")
    print(f"P_grid = {P_grid_final:.1f} MW")

solve_power_optimization()