import math

# --- Plan ---
# 1. Define given values and constants.
# 2. Convert units from mm to m.
# 3. Apply the formula for capacitance in F/m.
# 4. Convert the result to µF/km.
# 5. Print the full equation with numbers and the final result.

# --- Step 1: Define values ---
# Given parameters
r_w_mm = 11.25  # Wire radius in mm
m_mm = 17.32    # Distance from wire center to cable center in mm
R_mm = 32.32    # Inner radius of the screen in mm
e_r = 4.2       # Relative permittivity of the insulator

# Constants
e_0 = 8.854e-12 # Permittivity of free space in F/m

# --- Step 2: Convert units to meters ---
r_w = r_w_mm / 1000.0
m = m_mm / 1000.0
R = R_mm / 1000.0

# --- Step 3: Calculate capacitance in F/m ---
# Formula: C = (2 * pi * e_0 * e_r) / ln((R^2 - m^2) / (r_w * R))
try:
    log_argument = (R**2 - m**2) / (r_w * R)
    if log_argument <= 1:
        print("Error: Invalid physical dimensions. The argument of the logarithm must be greater than 1.")
    else:
        C_F_per_m = (2 * math.pi * e_0 * e_r) / math.log(log_argument)

        # --- Step 4: Convert to µF/km ---
        # To convert F/m to µF/km, we multiply by 10^9.
        C_uF_per_km = C_F_per_m * 1e9

        # --- Step 5: Print the equation and result ---
        print("The equation for capacitance (C) in µF/km is:")
        print(f"C = (2 * pi * {e_0} * {e_r} / ln(( {R}**2 - {m}**2 ) / ( {r_w} * {R} ))) * 10^9")
        print("\nResult:")
        print(f"C = {C_uF_per_km:.3f} µF/km")
        
        # Final answer for the system
        final_answer = f"<<<{C_uF_per_km:.3f}>>>"

except ValueError:
    print("Error during calculation. Please check the input values.")
except ZeroDivisionError:
    print("Error: Division by zero. Ensure r_w and R are not zero.")
