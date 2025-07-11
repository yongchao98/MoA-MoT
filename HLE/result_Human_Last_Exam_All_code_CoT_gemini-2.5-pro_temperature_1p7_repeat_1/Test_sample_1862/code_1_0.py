# --- Problem Parameters ---

# Case 1: Initial full sphere
diameter1_nm = 50.0  # Diameter in nanometers
time1_s = 20.0       # Lithiation time in seconds

# Case 2: Final quarter-sphere
# This particle is described as a quarter of a sphere with a 250 nm diameter.
# This means it was cut from an original sphere of 250 nm diameter.
diameter2_nm = 250.0 # Diameter of the original sphere in nanometers

# --- Analysis and Calculation ---

# Step 1: Calculate the radius for each case.
# The characteristic length for the lithiation process is the radius,
# which is the distance the phase-change front must travel.
radius1_nm = diameter1_nm / 2
radius2_nm = diameter2_nm / 2

# Step 2: Apply the physical model.
# The lithiation of LFP is modeled as a phase boundary moving at a constant
# velocity (v) from the surface to the center. Thus, time (t) is directly
# proportional to the radius (r), following the relationship t = r / v.
# From this, we can derive the ratio: time2 / time1 = radius2 / radius1.
# We can then solve for time2.
time2_s = time1_s * (radius2_nm / radius1_nm)

# Step 3: Display the results.
# The user wants to see the final equation with all numbers included.
print("Physical Model: Lithiation time (t) is proportional to particle radius (r).")
print("Formula: t_new = t_initial * (r_new / r_initial)")
print("\nCalculation:")
print(f"{time2_s:.1f} s = {time1_s:.1f} s * ({radius2_nm:.1f} nm / {radius1_nm:.1f} nm)")