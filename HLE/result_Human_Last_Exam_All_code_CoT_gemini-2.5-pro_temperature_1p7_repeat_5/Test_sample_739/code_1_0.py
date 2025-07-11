import math

# Plan: Calculate water film thickness (WFT) on the road surface.
# This will be done in two steps:
# 1. Calculate the flow rate per unit width (q) based on rainfall intensity and flow path length.
# 2. Use Manning's equation for shallow flow to calculate the water depth (d).

# --- Step 1: Define Parameters ---

# Rainfall Intensity (I) in mm/hr.
# Since IDF curves are not provided, we assume a standard design value for hydroplaning analysis,
# representing a high-intensity, short-duration storm.
I = 150.0  # mm/hr

# Flow Path Length (L) in meters.
# Water flows over 3 lanes, each 3.6 m wide.
num_lanes = 3
lane_width = 3.6  # m
L = num_lanes * lane_width  # m

# Pavement Cross-Slope (S) as a decimal.
# Given as 1.75%.
S = 1.75 / 100.0  # m/m

# Manning's Roughness Coefficient (n).
# For rough-textured asphalt, a typical value is 0.015.
n = 0.015

# --- Step 2: Perform Calculations ---

# 1. Calculate flow rate per unit width (q) in m^2/s.
# The formula is q = (I * L) / (3.6 * 10^6), which converts I from mm/hr to m/s.
q = (I * L) / (3.6 * 10**6)

# 2. Calculate water film thickness (d) in meters using the Izzard/Manning formula.
# The formula is: d = ( (q * n) / S^0.5 )^0.6
numerator = q * n
denominator = S**0.5
d_meters = (numerator / denominator)**0.6

# Convert the result from meters to millimeters.
d_mm = d_meters * 1000

# --- Step 3: Print the results ---

print("--- Calculation of Design Water Film Thickness ---")
print(f"Assumed and Given Parameters:")
print(f"  - Rainfall Intensity (I): {I} mm/hr (Assumed for hydroplaning analysis)")
print(f"  - Manning's Roughness (n): {n} (For rough-textured asphalt)")
print(f"  - Flow Path Length (L): {L:.1f} m ({num_lanes} lanes x {lane_width} m)")
print(f"  - Cross-Slope (S): {S} m/m")
print("\nFirst, we calculate the unit flow rate (q):")
print(f"q = (I * L) / (3.6 * 10^6)")
print(f"q = ({I} * {L:.1f}) / (3.6 * 10^6) = {q:.6f} m^2/s")

print("\nNext, we calculate the water film thickness (d) in mm:")
print("The formula used is: d_mm = [ (q * n) / (S**0.5) ]**0.6 * 1000")
# Output the final equation with all numbers substituted, as requested.
print("\nFinal Equation:")
print(f"d_mm = [ ({q:.6f} * {n}) / ({S}**0.5) ]**0.6 * 1000")
print(f"d_mm = [ ({numerator:.8f}) / ({denominator:.5f}) ]**0.6 * 1000")
print(f"d_mm = [ {(numerator / denominator):.8f} ]**0.6 * 1000")
print(f"\nResult:")
print(f"The design water film thickness is {d_mm:.2f} mm.")

# The final answer in the required format will be extracted from this result.
final_answer = round(d_mm, 2)