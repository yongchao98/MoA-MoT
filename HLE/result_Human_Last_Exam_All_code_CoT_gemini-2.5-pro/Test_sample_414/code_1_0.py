import numpy as np

def solve_radiochemistry_problem():
    """
    Calculates the time between sample irradiation and the first analysis
    based on the grow-in of Yttrium-90 from Strontium-90.
    """
    # --- Given data ---
    # Half-life of Yttrium-90 in days
    T_half_Y_hours = 64.1
    T_half_Y = T_half_Y_hours / 24.0

    # Activity measurements in kBq/mL
    A1 = 1.4
    A2 = 2.1

    # Time between measurements in days
    delta_t = 14.0

    # --- Calculations ---
    # Calculate the decay constant (lambda) for Y-90
    lambda_Y = np.log(2) / T_half_Y

    # We need to solve the equation: A2/A1 = (1 - exp(-l*(T+dt))) / (1 - exp(-l*T)) for T.
    # Let ratio = A2/A1, l = lambda_Y, dt = delta_t
    # ratio * (1 - exp(-l*T)) = 1 - exp(-l*T) * exp(-l*dt)
    # Let x = exp(-l*T)
    # ratio * (1 - x) = 1 - x * exp(-l*dt)
    # ratio - ratio*x = 1 - x * exp(-l*dt)
    # ratio - 1 = x * (ratio - exp(-l*dt))
    # x = (ratio - 1) / (ratio - exp(-l*dt))
    # exp(-l*T) = (ratio - 1) / (ratio - exp(-l*dt))
    # -l*T = log((ratio - 1) / (ratio - exp(-l*dt)))
    # T = -log((ratio - 1) / (ratio - exp(-l*dt))) / l

    ratio = A2 / A1
    
    # Calculate the term exp(-lambda * delta_t)
    exp_term = np.exp(-lambda_Y * delta_t)
    
    # Calculate the value of x = exp(-lambda * T)
    x = (ratio - 1) / (ratio - exp_term)
    
    # Calculate T
    T = -np.log(x) / lambda_Y

    # --- Output the results ---
    print("--- Problem Setup ---")
    print(f"Activity at T (A1): {A1} kBq/mL")
    print(f"Activity at T + 14 days (A2): {A2} kBq/mL")
    print(f"Y-90 Half-life: {T_half_Y:.4f} days")
    print(f"Y-90 Decay Constant (lambda): {lambda_Y:.4f} day^-1")
    print("\n--- Solving for T ---")
    print(f"The equation to solve is: {ratio:.1f} = (1 - exp(-{lambda_Y:.4f} * (T + {delta_t}))) / (1 - exp(-{lambda_Y:.4f} * T))")
    print(f"This simplifies to: exp(-{lambda_Y:.4f} * T) = ({ratio:.1f} - 1) / ({ratio:.1f} - exp(-{lambda_Y:.4f} * {delta_t}))")
    print(f"exp(-{lambda_Y:.4f} * T) = ({ratio - 1:.1f}) / ({ratio:.1f} - {exp_term:.4f}) = {x:.4f}")
    print(f"T = -ln({x:.4f}) / {lambda_Y:.4f}")
    print("\n--- Final Answer ---")
    print(f"The approximate time between irradiation and the first analysis is {T:.2f} days.")

solve_radiochemistry_problem()
<<<4.17>>>