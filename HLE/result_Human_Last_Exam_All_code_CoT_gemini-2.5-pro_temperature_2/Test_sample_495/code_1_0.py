import math

def solve_energy_loss():
    """
    Calculates the energy loss per centimetre for an alpha-particle based on Geiger's rule.
    """
    # --- Given parameters ---
    E0 = 8.5      # Initial energy in MeV
    R_total = 8.3 # Total range in cm
    x = 4.0       # Distance from the source in cm

    # --- Calculation Steps ---

    # Based on Geiger's rule, R = a * E^(3/2), we derive the formula for energy loss -dE/dx.
    # 1. Express E as a function of remaining range (R_total - x): E(x) = ((R_total - x) / a)^(2/3)
    # 2. Differentiate E(x) w.r.t x: dE/dx = -(2/3) * (1/a)^(2/3) * (R_total - x)^(-1/3)
    # 3. Energy loss is -dE/dx.

    # First, calculate the constant 'a' from the initial conditions: R_total = a * E0^(3/2)
    a = R_total / (E0 ** 1.5)

    # Then, calculate the energy loss per centimeter, -dE/dx, at distance x.
    # Breaking down the formula: -dE/dx = (2/3) * (1/a)^(2/3) * (R_total - x)^(-1/3)
    energy_loss_per_cm = (2/3) * (1/a)**(2/3) * (R_total - x)**(-1/3)

    # --- Final Output ---
    print(f"To find the energy loss per centimetre for an α-particle, we use Geiger's rule (R = a * E^(3/2)).")
    print(f"Given: Initial Energy E₀ = {E0} MeV, Total Range R = {R_total} cm, at distance x = {x} cm.")
    
    print("\nStep 1: Calculate the proportionality constant 'a'.")
    print(f"a = R / E₀^(3/2)")
    print(f"a = {R_total} / {E0}^(1.5) = {a:.4f}")

    print("\nStep 2: Calculate the energy loss per centimetre (-dE/dx) using the derived formula.")
    print(f"-dE/dx = (2/3) * (1/a)^(2/3) * (R - x)^(-1/3)")
    print(f"-dE/dx = (2/3) * (1/{a:.4f})^(2/3) * ({R_total} - {x})^(-1/3)")
    
    term1 = 2/3
    term2 = (1/a)**(2/3)
    term3 = (R_total - x)**(-1/3)
    
    print(f"-dE/dx = {term1:.4f} * {term2:.4f} * {term3:.4f}")

    print(f"\nThe final calculated energy loss per centimetre at {x} cm is {energy_loss_per_cm:.3f} MeV/cm.")

solve_energy_loss()