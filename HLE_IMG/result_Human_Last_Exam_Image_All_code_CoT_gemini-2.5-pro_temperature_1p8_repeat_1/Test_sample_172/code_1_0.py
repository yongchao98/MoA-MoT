import cmath

def solve_power_optimization():
    """
    Solves the hybrid power system optimization problem.
    """
    # Given values from the problem description
    P_wind_0 = 10.0  # MW (Initial wind power)
    P_wave_max = 5.0  # MW (Maximum wave power)
    P_grid_min = 12.0  # MW (Minimum required grid power)
    Pc_min = -2.0  # MW (Min compensation power)
    Pc_max = 3.0   # MW (Max compensation power)
    k = 0.01         # Loss coefficient
    V = 1.0          # p.u. (Bus voltage)

    # --- Step 1: Formulate the power equations ---
    # The total power delivered to the grid is P_grid = P_bus - P_loss.
    # The power loss is P_loss = k * I^2. With V=1, I=P_bus, so P_loss = k * P_bus^2.
    # The optimization constraint is P_bus - k*P_bus^2 >= P_grid_min.

    # --- Step 2: Solve for the optimal bus power (P_bus_opt) ---
    # To minimize loss (k*P_bus^2), we must minimize P_bus. We find the minimum
    # P_bus by solving the quadratic equation k*P_bus^2 - P_bus + P_grid_min = 0.
    a = k
    b = -1
    c = P_grid_min
    
    # Using the quadratic formula to find the roots
    discriminant = (b**2) - 4*(a*c)
    
    # Two potential solutions for P_bus
    sol1 = (-b - cmath.sqrt(discriminant)) / (2*a)
    sol2 = (-b + cmath.sqrt(discriminant)) / (2*a)
    
    # To minimize P_bus (and thus losses), we choose the smaller positive real root.
    P_bus_opt = min(sol1.real, sol2.real)

    # --- Step 3: Determine required wave power (P_wave_opt) ---
    # From the diagram, P_bus = P_wind_0 + P_wave.
    # So, P_wave_opt = P_bus_opt - P_wind_0
    P_wave_opt = P_bus_opt - P_wind_0

    # --- Step 4: Determine optimal compensation power (P_c_opt) ---
    # The total power and losses are independent of Pc. Assuming a secondary
    # objective to minimize internal control effort, we choose Pc_opt = 0.
    Pc_opt = 0.0

    # --- Step 5: Calculate final results ---
    P_loss_opt = k * P_bus_opt**2
    P_grid_final = P_bus_opt - P_loss_opt

    print("Formulation of the Total Power Delivered to the Grid:")
    print(f"The total power at the system bus is P_bus = P_wind^0 + P_wave = {P_wind_0:.2f} + P_wave.")
    print(f"The power loss in the cable is P_loss = k * P_bus^2 = {k:.2f} * P_bus^2.")
    print(f"The power delivered to the grid is P_grid = P_bus - P_loss.")
    print(f"The constraint is P_grid >= {P_grid_min:.2f} MW.\n")
    
    print("Optimization and Solution:")
    print(f"To minimize loss while meeting the demand, the target bus power is P_bus_opt = {P_bus_opt:.2f} MW.")
    print(f"This requires the wave sector to generate P_wave_opt = {P_wave_opt:.2f} MW.")
    print(f"The optimal compensation power is Pc_opt = {Pc_opt:.2f} MW.\n")

    print("Final Calculation for Power Delivered to the Grid:")
    equation_str = (
        f"P_grid = (P_wind^0 + P_wave_opt) - k * (P_wind^0 + P_wave_opt)^2\n"
        f"       = ({P_wind_0:.2f} + {P_wave_opt:.2f}) - {k:.2f} * ({P_wind_0:.2f} + {P_wave_opt:.2f})^2\n"
        f"       = ({P_bus_opt:.2f}) - {k:.2f} * ({P_bus_opt:.2f})^2\n"
        f"       = {P_bus_opt:.2f} - {P_loss_opt:.2f}\n"
        f"       = {P_grid_final:.2f} MW"
    )
    print(equation_str)

    print("\n--- Summary of Results ---")
    print(f"Optimal Compensation Power (Pc_opt): {Pc_opt:.2f} MW")
    print(f"Total Power Delivered to Grid: {P_grid_final:.2f} MW")

# Execute the function to solve the problem
solve_power_optimization()
>>>