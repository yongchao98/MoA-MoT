import math

def solve_power_optimization():
    """
    Calculates the optimal compensation power to minimize cable losses while meeting grid demand.
    """
    # System parameters
    P_wave_max = 5.0  # MW, maximum wave power
    P0_wind = 10.0    # MW, initial wind power
    P_grid_min = 12.0 # MW, minimum power to deliver to the grid
    Pc_min = -2.0     # MW, minimum compensation power
    Pc_max = 3.0      # MW, maximum compensation power
    k = 0.01          # Loss coefficient
    V = 1.0           # p.u., Bus voltage

    # We assume the analysis is for the case where the wave sector provides its maximum power.
    P_wave = P_wave_max

    # Step 1: Solve the grid power constraint to find the valid range for P_bus.
    # The constraint is P_grid >= P_grid_min, which is:
    # P_bus - k * P_bus^2 >= P_grid_min
    # This can be rewritten as a quadratic inequality:
    # k*P_bus^2 - P_bus + P_grid_min <= 0
    # Let's find the roots of the quadratic equation a*x^2 + b*x + c = 0 where x = P_bus.
    a = k
    b = -1
    c = P_grid_min
    
    # Using the quadratic formula: x = [-b +/- sqrt(b^2 - 4ac)] / 2a
    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        print("Error: No real solution for P_bus. The grid demand can never be met.")
        return

    sqrt_discriminant = math.sqrt(discriminant)
    P_bus_root1 = (-b - sqrt_discriminant) / (2 * a)
    P_bus_root2 = (-b + sqrt_discriminant) / (2 * a)

    # Since the parabola opens upwards (a > 0), the inequality holds between the roots.
    P_bus_lower_bound = P_bus_root1
    P_bus_upper_bound = P_bus_root2

    # Step 2: Determine the optimal P_c.
    # To minimize loss (k * P_bus^2), we must minimize P_bus.
    # The minimum P_bus that satisfies the grid constraint is P_bus_lower_bound.
    P_bus_opt = P_bus_lower_bound
    
    # Now find the required Pc to achieve this P_bus_opt.
    # P_bus = P0_wind + P_wave + Pc => Pc = P_bus - P0_wind - P_wave
    Pc_opt_candidate = P_bus_opt - P0_wind - P_wave

    # Check if this optimal Pc is within the allowed bounds [-2, 3].
    # If Pc_opt_candidate is smaller than Pc_min, we must use Pc_min.
    # If Pc_opt_candidate is larger than Pc_max, the problem setup is contradictory, but we'll cap it at Pc_max.
    # In our case, minimizing P_bus means minimizing Pc. We need the smallest Pc that is still feasible.
    # The feasible range for Pc is defined by both its own limits and the limits derived from P_bus bounds.
    
    Pc_from_Pbus_lower = P_bus_lower_bound - P0_wind - P_wave
    Pc_from_Pbus_upper = P_bus_upper_bound - P0_wind - P_wave
    
    # Combine constraints: max(Pc_min, Pc_from_Pbus_lower) <= Pc <= min(Pc_max, Pc_from_Pbus_upper)
    final_Pc_lower_bound = max(Pc_min, Pc_from_Pbus_lower)
    final_Pc_upper_bound = min(Pc_max, Pc_from_Pbus_upper)
    
    # To minimize loss, we choose the smallest possible Pc from the final valid range.
    Pc_opt = final_Pc_lower_bound

    # Step 3: Calculate the final power values with the optimal Pc.
    P_bus_final = P0_wind + P_wave + Pc_opt
    P_loss_final = k * (P_bus_final**2)
    P_grid_final = P_bus_final - P_loss_final

    # Print the results
    print(f"Given values:")
    print(f"  Initial Wind Power (P^0_wind): {P0_wind} MW")
    print(f"  Maximum Wave Power (P_wave): {P_wave} MW")
    print(f"  Minimum Grid Delivery (P_grid_min): {P_grid_min} MW")
    print(f"  Compensation Power (P_c) range: [{Pc_min}, {Pc_max}] MW")
    print(f"  Loss Coefficient (k): {k}\n")

    print(f"Optimization Results:")
    print(f"  To minimize losses, the system must generate the minimum power required to meet the 12 MW demand.")
    print(f"  The optimal compensation power (Pc_opt) is: {Pc_opt:.3f} MW\n")
    
    print(f"Final Power Calculation:")
    print(f"  Total power at the bus (P_bus) = {P0_wind} + {P_wave} + ({Pc_opt:.3f}) = {P_bus_final:.3f} MW")
    print(f"  Power loss in cable (P_loss) = {k} * ({P_bus_final:.3f})^2 = {P_loss_final:.3f} MW")
    print(f"  Total power delivered to grid (P_grid) = {P_bus_final:.3f} - {P_loss_final:.3f} = {P_grid_final:.3f} MW\n")

    print("Final equation for power delivered to the grid:")
    print(f"P_grid = (P^0_wind + P_wave + Pc_opt) - k * (P^0_wind + P_wave + Pc_opt)^2")
    print(f"{P_grid_final:.3f} = ({P0_wind} + {P_wave} + ({Pc_opt:.3f})) - {k} * ({P0_wind} + {P_wave} + ({Pc_opt:.3f}))^2")
    print(f"{P_grid_final:.3f} = {P_bus_final:.3f} - {k} * ({P_bus_final:.3f})^2")

solve_power_optimization()
<<< -1.055 >>>