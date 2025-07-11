import numpy as np

# --- 1. Define constants and initial values ---
P_wind_initial = 10  # MW
P_wave = 5           # MW
P_grid_min = 12      # MW
Pc_min = -2          # MW
Pc_max = 3           # MW
k = 0.01             # Loss coefficient in 1/MW (assuming P in MW, V=1)

# --- 2. Formulate and solve the optimization for P_bus ---
# Objective: Minimize P_loss = k * P_bus^2, which is equivalent to minimizing P_bus.
# Constraint: P_grid >= P_grid_min
# P_grid = P_bus - P_loss = P_bus - k * P_bus^2
# The constraint becomes: P_bus - k * P_bus^2 >= P_grid_min
# Rearranging into standard quadratic form (a*x^2 + b*x + c <= 0):
# k * P_bus^2 - P_bus + P_grid_min <= 0

# --- 3. Solve the quadratic inequality for P_bus ---
# Coefficients of the quadratic equation a*x^2 + b*x + c = 0
a = k
b = -1
c = P_grid_min

# Calculate the discriminant
delta = b**2 - 4 * a * c

if delta < 0:
    print("Error: The grid demand cannot be met as there are no real solutions for the required power.")
else:
    # Find the roots of the quadratic equation.
    # These roots define the boundaries of the feasible region for P_bus.
    P_bus_root1 = (-b - np.sqrt(delta)) / (2 * a)
    P_bus_root2 = (-b + np.sqrt(delta)) / (2 * a)

    # Since the parabola (a*x^2...) opens upwards (a>0), the inequality is satisfied between the roots.
    P_bus_feasible_min = P_bus_root1
    P_bus_feasible_max = P_bus_root2

    # To minimize losses, we must choose the minimum possible P_bus from the feasible range.
    P_bus_opt = P_bus_feasible_min

    # --- 4. Determine the optimal compensation power, Pc_opt ---
    # Based on the interpretation P_bus = (P_wind_initial + P_wave) - Pc
    # P_bus_opt = (P_wind_initial + P_wave) - Pc_opt
    Pc_opt = (P_wind_initial + P_wave) - P_bus_opt

    # --- 5. Verify Pc_opt is within its bounds and calculate final values ---
    if Pc_min <= Pc_opt <= Pc_max:
        P_loss_opt = k * P_bus_opt**2
        P_grid_opt = P_bus_opt - P_loss_opt

        print("Formulation:")
        print(f"Total Power at Bus: P_bus = P_wind^0 + P_wave - P_c = ({P_wind_initial} + {P_wave}) - P_c")
        print(f"Power Loss: P_loss = k * P_bus^2 = {k} * P_bus^2")
        print(f"Power to Grid: P_grid = P_bus - P_loss")
        print("-" * 30)

        print("Optimization Results:")
        print(f"To deliver at least {P_grid_min} MW, the feasible range for bus power is [{P_bus_feasible_min:.3f}, {P_bus_feasible_max:.3f}] MW.")
        print(f"To minimize losses, the optimal bus power is chosen as the minimum of this range.")
        print(f"Optimal Power at Bus (P_bus_opt): {P_bus_opt:.3f} MW")
        print("-" * 30)

        print("Final Calculated Values:")
        print(f"Optimal Compensation Power (Pc_opt): {Pc_opt:.3f} MW")
        print(f"Resulting Power Loss (P_loss_opt): {P_loss_opt:.3f} MW")
        print(f"Total Power Delivered to Grid (P_grid_opt): {P_grid_opt:.3f} MW")
        print("-" * 30)

        print("Final Power Delivery Equation:")
        print("P_grid = P_bus - P_loss")
        print(f"{P_grid_opt:.3f} MW = {P_bus_opt:.3f} MW - {P_loss_opt:.3f} MW")

    else:
        print(f"Error: The calculated optimal Pc ({Pc_opt:.3f} MW) is outside the allowed bounds [{Pc_min}, {Pc_max}] MW.")

<<<1.056>>>