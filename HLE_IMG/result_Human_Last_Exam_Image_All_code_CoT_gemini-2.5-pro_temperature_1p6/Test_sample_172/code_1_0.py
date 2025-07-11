import math

def solve_power_optimization():
    """
    Solves for the optimal compensation power to minimize cable losses
    while meeting grid demand.
    """
    # Given parameters
    P0_wind = 10.0  # MW
    P_wave = 5.0    # MW (operating at max)
    P_grid_min = 12.0 # MW
    Pc_min = -2.0   # MW
    Pc_max = 3.0    # MW
    k = 0.01        # Loss coefficient (1/MW for consistency)
    V_bus = 1.0     # p.u.

    # Step 1: Formulate the optimization problem in terms of P_bus
    # P_bus = P0_wind + P_wave - P_c = 15 - P_c
    # P_loss = k * P_bus^2
    # P_grid = P_bus - P_loss >= P_grid_min
    # We need to solve for P_bus in: P_bus - k * P_bus^2 >= P_grid_min
    # k*P_bus^2 - P_bus + P_grid_min <= 0
    # 0.01 * P_bus^2 - P_bus + 12 <= 0

    # Step 2: Solve the quadratic equation 0.01*x^2 - x + 12 = 0 to find the boundaries for P_bus
    a = k
    b = -1.0
    c = P_grid_min
    
    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        print("No real solution for P_bus exists. The grid demand can never be met.")
        return

    sqrt_discriminant = math.sqrt(discriminant)
    
    # The valid region for P_bus is between the two roots of the quadratic equation
    P_bus_sol1 = (-b - sqrt_discriminant) / (2 * a)
    P_bus_sol2 = (-b + sqrt_discriminant) / (2 * a)
    
    P_bus_min_req = min(P_bus_sol1, P_bus_sol2)
    P_bus_max_req = max(P_bus_sol1, P_bus_sol2)
    
    # Step 3: Determine the optimal P_bus
    # To minimize loss, we need to minimize P_bus.
    # The smallest P_bus that satisfies the grid constraint is P_bus_min_req.
    P_bus_opt = P_bus_min_req

    # Step 4: Calculate the optimal compensation power Pc_opt
    # P_bus_opt = 15 - Pc_opt  => Pc_opt = 15 - P_bus_opt
    P_initial_total = P0_wind + P_wave
    Pc_opt = P_initial_total - P_bus_opt
    
    # Check if Pc_opt is within its allowed range [-2, 3]
    if not (Pc_min <= Pc_opt <= Pc_max):
        # This part of the code would handle cases where the optimal Pc is outside its bounds.
        # In this specific problem, the calculated Pc_opt falls within the bounds.
        # If it were outside, we'd choose the bound (Pc_min or Pc_max) that results in the
        # lowest P_bus while still satisfying P_grid >= 12.
        # For this problem, our calculated Pc_opt (approx 1.056) is valid.
        pass

    # Step 5: Calculate the final power values using the optimal P_c
    P_loss_opt = k * P_bus_opt**2
    P_grid_opt = P_bus_opt - P_loss_opt

    # Print the results
    print("--- System Power Optimization ---")
    print(f"Initial Wind Power (P^0_wind): {P0_wind} MW")
    print(f"Wave Power (P_wave): {P_wave} MW")
    print(f"Minimum Required Grid Power (P_grid_min): {P_grid_min} MW")
    
    print("\n--- Optimal Power Configuration ---")
    print(f"Optimal Compensation Power (Pc_opt): {Pc_opt:.3f} MW")
    print(f"This value is found by selecting the minimum possible bus power that satisfies the grid demand.")
    
    print("\n--- Resulting Power Flow ---")
    print(f"Total Power at System Bus (P_bus): {P_bus_opt:.3f} MW")
    print(f"Undersea Cable Power Loss (P_loss): {P_loss_opt:.3f} MW")
    print(f"Total Power Delivered to Grid (P_grid): {P_grid_opt:.3f} MW")

    print("\nFinal Power Delivery Equation:")
    print(f"P_grid = P_bus - P_loss")
    print(f"{P_grid_opt:.3f} MW = {P_bus_opt:.3f} MW - {P_loss_opt:.3f} MW")


solve_power_optimization()
<<<1.056>>>