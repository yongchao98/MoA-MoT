import math

def solve_power_system():
    """
    Solves for the optimal compensation power and grid power delivery
    in a hybrid offshore power generation system.
    """
    # System parameters
    P_wave = 5.0          # Maximum wave power in MW
    P0_wind = 10.0        # Initial wind power in MW
    P_grid_demand = 12.0  # Minimum power to be delivered to the grid in MW
    Pc_min = -2.0         # Minimum compensation power in MW
    Pc_max = 3.0          # Maximum compensation power in MW
    k = 0.01              # Loss coefficient
    V = 1.0               # Bus voltage in p.u.

    # --- Step 1 & 2: Formulate Equations and Objective ---
    print("Step 1: Formulating the Power Equations and Objective")
    print("The total power injected into the bus (P_bus) is modeled as:")
    print(f"P_bus = P^0_wind + P_wave + P_c = {P0_wind} + {P_wave} + P_c")
    print("The power loss in the cable is P_loss = k * P_bus^2.")
    print("The power delivered to the grid is P_grid = P_bus - P_loss.")
    print("The objective is to minimize P_loss, which means minimizing P_bus while meeting all constraints.\n")

    # --- Step 3: Analyze Grid Demand Constraint ---
    print("Step 2: Solving the Grid Demand Constraint")
    print(f"The grid demand requires: P_bus - P_loss >= {P_grid_demand} MW")
    print(f"Substituting for P_loss: P_bus - {k} * P_bus^2 >= {P_grid_demand}")
    print(f"This can be written as a quadratic inequality: {k} * P_bus^2 - P_bus + {P_grid_demand} <= 0")
    
    # Solve the quadratic equation k*x^2 - x + P_grid_demand = 0 to find the boundaries
    a = k
    b = -1.0
    c = P_grid_demand
    discriminant = b**2 - 4 * a * c
    
    if discriminant < 0:
        print("Error: The demand cannot be met as the quadratic has no real roots.")
        return

    # The valid range for P_bus is between the two roots.
    # To minimize P_bus, we choose the smaller root.
    P_bus_min_req = (-b - math.sqrt(discriminant)) / (2 * a)
    print(f"Solving this gives the minimum required bus power to meet demand: P_bus_req = {P_bus_min_req:.4f} MW\n")

    # --- Step 4: Determine Optimal Compensation Power ---
    print("Step 3: Determining the Optimal Compensation Power (P_c^opt)")
    P_base = P0_wind + P_wave
    Pc_req = P_bus_min_req - P_base
    print(f"To achieve this minimum bus power, the required P_c is:")
    print(f"P_c_req = P_bus_req - (P^0_wind + P_wave) = {P_bus_min_req:.4f} - ({P0_wind} + {P_wave}) = {Pc_req:.4f} MW")
    
    # The optimal P_c must satisfy both its own bounds and the demand requirement.
    # We need P_c >= Pc_req and P_c >= Pc_min. To minimize P_c, we take the max of the lower bounds.
    Pc_opt = max(Pc_req, Pc_min)
    
    print(f"P_c is also bounded between {Pc_min} MW and {Pc_max} MW.")
    print(f"To minimize P_bus, we need the smallest valid P_c, which is max({Pc_req:.4f}, {Pc_min}).")
    
    if Pc_opt > Pc_max:
        print(f"Error: Optimal P_c ({Pc_opt:.4f} MW) is not achievable as it exceeds the maximum limit of {Pc_max} MW.")
        return
        
    print(f"Therefore, the optimal compensation power is P_c^opt = {Pc_opt:.4f} MW.\n")
    
    # --- Step 5: Calculate Final Power Delivery ---
    print("Step 4: Calculating Final Power Values and Grid Delivery")
    P_bus_opt = P_base + Pc_opt
    P_loss_opt = k * P_bus_opt**2
    P_grid_final = P_bus_opt - P_loss_opt
    
    print("Using P_c^opt, the final equation for power delivered to the grid is:")
    print(f"P_grid = (P^0_wind + P_wave + P_c^opt) - k * (P^0_wind + P_wave + P_c^opt)^2")
    print(f"P_grid = ({P0_wind} + {P_wave} + ({Pc_opt:.4f})) - {k} * ({P0_wind} + {P_wave} + ({Pc_opt:.4f}))^2")
    print(f"P_grid = ({P_bus_opt:.4f}) - {k} * ({P_bus_opt:.4f})^2")
    print(f"P_grid = {P_bus_opt:.4f} - {P_loss_opt:.4f}")
    print(f"\nFinal Result:")
    print(f"The total power delivered to the grid is {P_grid_final:.4f} MW.")
    print(f"The optimal compensation power required is {Pc_opt:.4f} MW.")

solve_power_system()
<<< -1.055 >>>