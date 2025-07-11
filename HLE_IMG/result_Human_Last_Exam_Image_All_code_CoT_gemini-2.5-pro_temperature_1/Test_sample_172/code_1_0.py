import math

def solve_power_optimization():
    """
    Solves the hybrid power system optimization problem.
    """
    # --- Given Parameters ---
    P_wave = 5.0  # MW
    P_wind_0 = 10.0  # MW
    P_grid_min = 12.0  # MW
    Pc_min = -2.0  # MW
    Pc_max = 3.0   # MW
    k = 0.01  # Loss coefficient
    V = 1.0  # p.u. bus voltage

    print("--- Problem Formulation ---")
    print("1. Total Power at the System Bus (P_bus):")
    print(f"P_bus = P_wind^0 + P_wave + Pc = {P_wind_0} + {P_wave} + Pc = {P_wind_0 + P_wave} + Pc")
    
    print("\n2. Power Loss in the Undersea Cable (P_loss):")
    print(f"P_loss = k * (P_bus / V)^2 = {k} * P_bus^2")

    print("\n3. Power Delivered to the Onshore Grid (P_grid):")
    print("P_grid = P_bus - P_loss")
    print(f"P_grid = P_bus - {k} * P_bus^2\n")

    # --- Solving the Optimization Problem ---
    # Objective: Minimize P_loss, which means minimizing P_bus.
    # Constraint 1: P_grid >= P_grid_min
    # P_bus - k * P_bus^2 >= P_grid_min  =>  k*P_bus^2 - P_bus + P_grid_min <= 0
    # We solve the quadratic equation k*x^2 - x + P_grid_min = 0 to find the boundaries.
    a = k
    b = -1
    c = P_grid_min
    
    # Calculate the discriminant
    discriminant = b**2 - 4 * a * c
    
    if discriminant < 0:
        print("Error: The grid demand cannot be met with the given parameters.")
        return

    # Calculate the roots for P_bus
    sqrt_discriminant = math.sqrt(discriminant)
    p_bus_root1 = (-b - sqrt_discriminant) / (2 * a)
    p_bus_root2 = (-b + sqrt_discriminant) / (2 * a)
    
    # The feasible range for P_bus from the grid demand constraint
    p_bus_feasible_1_min = min(p_bus_root1, p_bus_root2)
    p_bus_feasible_1_max = max(p_bus_root1, p_bus_root2)

    # Constraint 2: Pc_min <= Pc <= Pc_max
    # This defines another feasible range for P_bus, since Pc = P_bus - (P_wind_0 + P_wave)
    p_bus_feasible_2_min = Pc_min + P_wind_0 + P_wave
    p_bus_feasible_2_max = Pc_max + P_wind_0 + P_wave

    # The final feasible range for P_bus is the intersection of the two ranges
    final_feasible_p_bus_min = max(p_bus_feasible_1_min, p_bus_feasible_2_min)
    final_feasible_p_bus_max = min(p_bus_feasible_1_max, p_bus_feasible_2_max)

    if final_feasible_p_bus_min > final_feasible_p_bus_max:
        print("Error: No feasible solution exists that satisfies all constraints.")
        return

    # To minimize loss, we choose the smallest possible P_bus in the feasible range
    P_bus_opt = final_feasible_p_bus_min

    # --- Calculate Optimal Values ---
    P_c_opt = P_bus_opt - (P_wind_0 + P_wave)
    P_loss_opt = k * P_bus_opt**2
    P_grid_opt = P_bus_opt - P_loss_opt
    
    print("--- Optimal Solution ---")
    print(f"The optimal compensation power required to minimize cable loss is:")
    print(f"Pc_opt = {P_c_opt:.3f} MW\n")

    print("The final power delivered to the grid is calculated with the following equation:")
    print(f"P_grid = (P_wind^0 + P_wave + Pc_opt) - k * (P_wind^0 + P_wave + Pc_opt)^2")
    print(f"P_grid = ({P_wind_0} + {P_wave} + {P_c_opt:.3f}) - {k} * ({P_wind_0} + {P_wave} + {P_c_opt:.3f})^2")
    print(f"P_grid = ({P_bus_opt:.3f}) - {k} * ({P_bus_opt:.3f})^2")
    print(f"P_grid = {P_bus_opt:.3f} - {P_loss_opt:.3f}")
    print(f"P_grid = {P_grid_opt:.3f} MW")

solve_power_optimization()
<<<P_grid = 13.945 - 1.945
P_grid = 12.000 MW>>>