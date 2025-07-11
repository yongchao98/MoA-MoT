import numpy as np

def solve_power_optimization():
    """
    Solves the hybrid power system optimization problem.
    """
    # --- Step 1: Define System Parameters ---
    P_wind_initial = 10  # MW
    P_wave_max = 5       # MW, assuming wave sector operates at maximum power for this scenario
    P_grid_min = 12      # MW
    Pc_min = -2          # MW
    Pc_max = 3           # MW
    k = 0.01             # Loss coefficient (p.u. resistance on a 1 MW or matched base)

    # For this problem, we assume the non-controllable wave power is at its maximum value
    # to find a specific numerical solution for the optimal compensation power.
    P_wave = P_wave_max

    # --- Step 2: Formulate Power Equations ---
    # The total power at the system bus (P_bus) is modeled as the sum of the initial
    # wind power, the wave power, and the controllable compensation power P_c.
    # The power loss in the cable is P_loss = k * P_bus^2.
    # The power delivered to the grid is P_grid = P_bus - P_loss.

    print("Step 1: Formulating the Power Equations")
    print("-" * 50)
    print("Based on the problem description, the key equations are:")
    print(f"Total Bus Power: P_bus = P^0_wind + P_wave + P_c")
    print(f"Cable Loss: P_loss = k * P_bus^2 = {k} * P_bus^2")
    print(f"Power to Grid: P_grid = P_bus - P_loss")
    print("\nUsing the given values (with P_wave at its maximum):")
    print(f"P_bus = {P_wind_initial} MW + {P_wave} MW + P_c = {P_wind_initial + P_wave} + P_c")
    print("-" * 50)


    # --- Step 3: Define and Solve the Optimization Problem ---
    # Objective: Minimize P_loss, which is equivalent to minimizing P_bus.
    # Subject to constraints:
    # 1. P_grid >= P_grid_min  =>  P_bus - k*P_bus^2 >= P_grid_min
    # 2. Pc_min <= P_c <= Pc_max

    print("\nStep 2: Solving the Optimization Problem")
    print("-" * 50)
    print("Objective: Minimize cable losses (P_loss) by minimizing total bus power (P_bus).")
    print(f"Constraint 1: Grid power must be at least {P_grid_min} MW.")
    print(f"Constraint 2: Compensation power P_c must be in [{Pc_min}, {Pc_max}] MW.")

    # Solve Constraint 1: P_bus - k*P_bus^2 >= P_grid_min
    # Rearranging gives a quadratic inequality: k*P_bus^2 - P_bus + P_grid_min <= 0
    coeffs = [k, -1, P_grid_min]
    roots = np.roots(coeffs)
    P_bus_min_req = min(roots)
    P_bus_max_req = max(roots)
    
    print("\nFrom Constraint 1 (Grid Demand):")
    print(f"Solving P_bus - {k}*P_bus^2 >= {P_grid_min} requires P_bus to be between the roots of the quadratic equation.")
    print(f"Required P_bus range: [{P_bus_min_req:.4f}, {P_bus_max_req:.4f}] MW")

    # Analyze Constraint 2: Pc_min <= P_c <= Pc_max
    # This defines a range for P_bus = P_wind_initial + P_wave + P_c
    P_bus_from_Pc_min = P_wind_initial + P_wave + Pc_min
    P_bus_from_Pc_max = P_wind_initial + P_wave + Pc_max

    print("\nFrom Constraint 2 (Compensation Power Limits):")
    print(f"The range for P_c translates to a P_bus range of [{P_bus_from_Pc_min}, {P_bus_from_Pc_max}] MW.")

    # Combine constraints to find the overall feasible range for P_bus
    feasible_P_bus_min = max(P_bus_min_req, P_bus_from_Pc_min)
    feasible_P_bus_max = min(P_bus_max_req, P_bus_from_Pc_max)
    
    print("\nThe overall feasible range for P_bus is the intersection of these two ranges:")
    print(f"Feasible P_bus range: [{feasible_P_bus_min:.4f}, {feasible_P_bus_max:.4f}] MW")

    # To minimize loss, we choose the minimum feasible P_bus.
    P_bus_opt = feasible_P_bus_min
    print(f"\nTo minimize loss, the optimal bus power is the minimum feasible value: {P_bus_opt:.3f} MW.")
    print("-" * 50)

    # --- Step 4: Calculate Final Optimal Values ---
    # Optimal compensation power P_c_opt
    P_c_opt = P_bus_opt - P_wind_initial - P_wave

    # Optimal loss and grid power
    P_loss_opt = k * P_bus_opt**2
    P_grid_final = P_bus_opt - P_loss_opt

    print("\nStep 3: Final Results")
    print("-" * 50)
    print(f"Optimal Compensation Power (P_c^opt):")
    print(f"P_c_opt = P_bus_opt - P^0_wind - P_wave = {P_bus_opt:.3f} - {P_wind_initial} - {P_wave} = {P_c_opt:.3f} MW")
    
    print(f"\nTotal Power Delivered to Grid (P_grid):")
    print("The final equation with optimal values inserted is:")
    print(f"P_grid = (P^0_wind + P_wave + P_c^opt) - k * (P^0_wind + P_wave + P_c^opt)^2")
    print(f"P_grid = ({P_wind_initial} + {P_wave} + ({P_c_opt:.3f})) - {k} * ({P_wind_initial} + {P_wave} + ({P_c_opt:.3f}))^2")
    print(f"P_grid = {P_bus_opt:.3f} - {k} * ({P_bus_opt:.3f})^2")
    print(f"P_grid = {P_bus_opt:.3f} - {P_loss_opt:.3f} = {P_grid_final:.3f} MW")
    print("-" * 50)

if __name__ == '__main__':
    solve_power_optimization()
    P_c_opt_val = 13.9445 - 10 - 5
    # The primary value asked for is the optimal compensation power
    print(f"\n<<<P_c^opt = {P_c_opt_val:.3f} MW>>>")