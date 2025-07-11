import math

def solve_power_system_optimization():
    """
    Solves the hybrid power system optimization problem.
    """
    # Step 1: Define given constants and constraints
    P0_wind = 10  # MW, initial combined wind power
    P_wave_max = 5  # MW, maximum wave power
    k = 0.01  # Loss coefficient
    V = 1.0  # p.u., bus voltage
    P_grid_min = 12  # MW, minimum power to be delivered to the grid
    Pc_min = -2  # MW, minimum compensation power
    Pc_max = 3  # MW, maximum compensation power

    # Assumption: The wave sector operates at its maximum power for this calculation.
    P_wave = P_wave_max

    print("Step-by-step Solution:")
    print("----------------------\n")

    # Step 2: Formulate the total power delivered to the grid
    print("1. Formulation of Power Delivered to the Grid (P_grid):")
    print("Based on the problem statement, we assume the total power at the bus is managed by Pc.")
    print("P_bus = P^0_wind + P_wave + Pc")
    print(f"P_bus = {P0_wind} MW + P_wave + Pc\n")
    print("Power loss is given by P_loss = k * (P_bus / V)^2.")
    print(f"P_loss = {k} * (P_bus / {V})^2\n")
    print("The power delivered to the grid is P_grid = P_bus - P_loss.")
    print("P_grid = (P^0_wind + P_wave + Pc) - k * ((P^0_wind + P_wave + Pc) / V)^2\n")

    # Step 3: Define the optimization problem
    print("2. Optimization Goal:")
    print("Minimize cable losses (P_loss) by minimizing the total power at the bus (P_bus),")
    print(f"subject to the constraint that P_grid >= {P_grid_min} MW.\n")

    # Step 4: Solve for the minimum required P_bus
    print("3. Solving for the optimal P_bus:")
    print(f"The constraint is: P_bus - {k} * P_bus^2 >= {P_grid_min}")
    print(f"Rearranging gives the quadratic inequality: {k}*P_bus^2 - P_bus + {P_grid_min} <= 0")
    
    # Solve the quadratic equation 0.01*P_bus^2 - 1*P_bus + 12 = 0
    a = k
    b = -1
    c = P_grid_min
    discriminant = b**2 - 4*a*c
    
    if discriminant < 0:
        print("The grid demand cannot be met.")
        return

    sqrt_discriminant = math.sqrt(discriminant)
    p_bus_root1 = (-b - sqrt_discriminant) / (2*a)
    p_bus_root2 = (-b + sqrt_discriminant) / (2*a)

    # The valid range for P_bus is between the two roots.
    # To minimize loss, we need the minimum possible P_bus.
    P_bus_opt = p_bus_root1
    print(f"The roots of the quadratic equation are {p_bus_root1:.3f} MW and {p_bus_root2:.3f} MW.")
    print(f"To meet the grid demand while minimizing losses, the optimal bus power is P_bus_opt = {P_bus_opt:.3f} MW.\n")

    # Step 5: Calculate the optimal compensation power Pc_opt
    print("4. Calculating the Optimal Compensation Power (Pc_opt):")
    print(f"Using the assumption that P_wave = {P_wave} MW (its maximum).")
    # P_bus_opt = P0_wind + P_wave + Pc_opt
    Pc_opt = P_bus_opt - P0_wind - P_wave
    print(f"Pc_opt = P_bus_opt - P^0_wind - P_wave = {P_bus_opt:.3f} - {P0_wind} - {P_wave} = {Pc_opt:.3f} MW.")

    # Step 6: Verify Pc_opt is within its bounds
    if not (Pc_min <= Pc_opt <= Pc_max):
        print(f"Warning: The calculated Pc_opt ({Pc_opt:.3f} MW) is outside the allowed range [{Pc_min}, {Pc_max}] MW.")
    else:
        print(f"The calculated Pc_opt ({Pc_opt:.3f} MW) is within the allowed range [{Pc_min}, {Pc_max}] MW.\n")
    
    # Step 7: Compute final power values
    P_loss_final = k * (P_bus_opt / V)**2
    P_grid_final = P_bus_opt - P_loss_final
    
    print("5. Final Results:")
    print("-----------------")
    print(f"Optimal Compensation Power (Pc_opt): {Pc_opt:.3f} MW")
    print(f"Total Power at Bus (P_bus): {P_bus_opt:.3f} MW")
    print(f"Undersea Cable Power Loss (P_loss): {P_loss_final:.3f} MW")
    print(f"Total Power Delivered to Grid (P_grid): {P_grid_final:.3f} MW\n")
    
    print("Final Power Balance Equation:")
    print("P_grid = P_bus - P_loss")
    print(f"{P_grid_final:.3f} MW = {P_bus_opt:.3f} MW - {P_loss_final:.3f} MW")

solve_power_system_optimization()
<<<
-1.056
>>>