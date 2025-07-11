import math

def analyze_titan_feasibility():
    """
    Analyzes the feasibility of calculating the rocket force problem on the Titan computer.
    """
    print("Starting analysis for the Titan Superconducting Computer...\n")

    # --- Step 1: Define the physics formula ---
    print("Step 1: Establishing the governing physical equation.")
    print("The horizontal distance `x` traveled under constant force `F` is x = (1/2) * a * t^2.")
    print("The horizontal acceleration `a` is F / m.")
    print("The time of flight `t` from height `y0` is t = sqrt(2 * y0 / g).")
    print("Combining these, we get: x = (1/2) * (F/m) * (2 * y0 / g)")
    print("Solving for the force F, the final equation is:")
    print("F = (x_final * m * g) / y0\n")
    print("where mass m = density * Volume = rho * (4/3) * pi * r^3.\n")


    # --- Step 2: Define problem constants in SI units ---
    print("Step 2: Listing the physical constants in standard SI units (meters, kg, seconds).")
    # Given values
    radius_cm = 0.5
    density_kg_per_cm3 = 0.9
    height_m = 10.0
    distance_m = 20.0 # We target the middle of the 19m-21m range
    gravity_ms2 = 9.8
    pi = math.pi

    # Convert to consistent SI units
    # r = 0.5 cm -> 0.005 m
    r_m = radius_cm / 100.0
    # rho = 0.9 kg/cm^3 -> 900,000 kg/m^3
    rho_kg_per_m3 = density_kg_per_cm3 * (100**3)

    constants = {
        "y0 (initial height)": height_m,
        "x_final (target distance)": distance_m,
        "g (gravity)": gravity_ms2,
        "pi": pi,
        "r (rock radius)": r_m,
        "rho (rock density)": rho_kg_per_m3,
    }
    for name, val in constants.items():
        print(f"- {name}: {val}")
    print("")

    # --- Step 3: Check if constants can be represented in Titan's 5-bit system ---
    print("Step 3: Checking if these constants can be represented as Titan 5-bit fractions.")
    print("A number must be represented as a/b, where a and b are integers from 0 to 31.")
    
    MAX_INT = 31
    MIN_REPRESENTABLE_POSITIVE = 1.0 / MAX_INT
    MAX_REPRESENTABLE = float(MAX_INT)
    
    print(f"This means any representable number must be within the approximate range [{MIN_REPRESENTABLE_POSITIVE:.4f}, {MAX_REPRESENTABLE:.1f}].\n")
    
    is_possible = True
    
    for name, value in constants.items():
        if value > MAX_REPRESENTABLE:
            print(f"Error: Constant '{name}' with value {value} is too large to be represented.")
            print(f"       It is greater than the maximum possible value of {MAX_REPRESENTABLE}.\n")
            is_possible = False
        elif 0 < value < MIN_REPRESENTABLE_POSITIVE:
            print(f"Error: Constant '{name}' with value {value} is too small to be represented.")
            print(f"       It is smaller than the minimum positive value of {MIN_REPRESENTABLE_POSITIVE:.4f}.\n")
            is_possible = False
        else:
            # Note: This is a simplified check. Even if within range, a suitable a/b might not exist.
            # But for this problem, the extreme values are the issue.
            print(f"OK: Constant '{name}' with value {value} is potentially representable.")
            if name == "g (gravity)":
                print("    (e.g., as 29/3 = 9.67)")
            if name == "pi":
                print("    (e.g., as 22/7 = 3.14)")

    # --- Step 4: Final Conclusion ---
    print("\n--- Final Analysis ---")
    if not is_possible:
        print("The calculation is NOT POSSIBLE.")
        print("As shown above, the rock's radius (r = 0.005 m) and density (rho = 900000.0 kg/m^3)")
        print("are far outside the numerical range of the Titan computer's 5-bit fractional system.")
        print("Since these fundamental inputs cannot be represented, the calculation cannot begin.")
    else:
        # This case is not reached in this problem
        print("The calculation is POSSIBLE.")

if __name__ == '__main__':
    analyze_titan_feasibility()
