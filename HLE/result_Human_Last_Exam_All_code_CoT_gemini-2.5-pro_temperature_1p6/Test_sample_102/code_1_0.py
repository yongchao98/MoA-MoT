import math

# Step 1: Define the given equivalent values
ni_eq = 29.0
cr_eq = 39.0

# Step 2: Define reference points from the WRC-1992 diagram for Cr_eq = 39
# These are the Nickel Equivalent values for constant Ferrite Number (FN) lines.
# FN is approximately equal to ferrite percentage.
fn_low = 60.0
ni_eq_for_fn_low = 29.6

fn_high = 70.0
ni_eq_for_fn_high = 28.4

# Step 3: Perform linear interpolation to find the ferrite level
# We calculate the proportion of the distance our Ni_eq point is along the line
# connecting the two reference points.
# Note: As Ni_eq increases, ferrite level decreases.
proportion = (ni_eq_for_fn_low - ni_eq) / (ni_eq_for_fn_low - ni_eq_for_fn_high)
ferrite_level = fn_low + proportion * (fn_high - fn_low)

# Step 4: Round the result to the nearest 10 (rounding .5 cases up)
rounded_ferrite_level = int((ferrite_level / 10) + 0.5) * 10

# Step 5: Print the explanation, calculation, and final answer
print(f"Given Nickel Equivalent (Ni_eq): {ni_eq}")
print(f"Given Chromium Equivalent (Cr_eq): {cr_eq}\n")
print("Using linear interpolation based on the WRC-1992 diagram at Cr_eq = 39:")
print(f"A Ferrite Number of {int(fn_low)} corresponds to Ni_eq of ~{ni_eq_for_fn_low}")
print(f"A Ferrite Number of {int(fn_high)} corresponds to Ni_eq of ~{ni_eq_for_fn_high}\n")

print("Calculating the unrounded ferrite level:")
print(f"Equation: {fn_low} + (({ni_eq_for_fn_low} - {ni_eq}) / ({ni_eq_for_fn_low} - {ni_eq_for_fn_high})) * ({fn_high} - {fn_low})")
print(f"Result: {fn_low} + {proportion:.2f} * {fn_high - fn_low} = {ferrite_level:.1f}\n")

print(f"Rounding the result ({ferrite_level:.1f}) to the nearest 10:")
print(f"Final Approximate Ferrite Level = {rounded_ferrite_level}")