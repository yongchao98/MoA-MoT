import math

# User-defined inputs
# n: represents the number of 100-liter water containers available at the origin.
# m: the distance to the destination in kilometers.
n = 3
m = 70

# --- Start of the solution logic ---

# Step 1: Find 'k', the number of trips required for the final segment of the journey.
# A phase 'j' is defined by requiring 'j' forward trips and 'j-1' return trips.
# The water consumption rate in this phase is (2*j - 1) L/km. A phase 'j' ends after
# traveling a distance of 100 / (2*j - 1) km, at which point one 100L load has been
# fully consumed, and the problem transitions to a phase with 'j-1' loads.

distance_covered_before_segment = 0.0
k = 0 # k will be the number of trips in the final segment.

# We loop from j=n down to 2 to see if the journey completes the segment.
for j in range(n, 1, -1):
    # This is the maximum distance that can be traveled in the phase with 'j' trips.
    max_segment_length = 100.0 / (2 * j - 1)
    
    # If the total distance 'm' falls within this segment.
    if m <= distance_covered_before_segment + max_segment_length:
        k = j
        break
    
    # If 'm' is greater, this segment is completed. Update the distance covered
    # and move to the next phase (with j-1 trips).
    distance_covered_before_segment += max_segment_length

# If the loop completes without finding k, it means 'm' is in the last phase
# where there is only one load left, requiring a single trip.
if k == 0:
    k = 1

# Step 2: Calculate the components for the final equation and construct the equation string.
# The general formula for water left is:
# Water Left = k*100 - (2*k - 1) * (m - distance_covered_before_segment)
# where distance_covered_before_segment = sum_{i=k+1 to n} (100 / (2*i - 1))

# Calculate the numerical value of the summation part.
summation_val = 0.0
for i in range(k + 1, n + 1):
    summation_val += 100.0 / (2 * i - 1)

# Calculate the final numerical answer for water left.
water_left = k * 100.0 - (2 * k - 1) * (m - summation_val)

# Construct the string for the summation part of the equation.
summation_str_parts = []
if k < n:
    for i in range(k + 1, n + 1):
        summation_str_parts.append(f"100/(2*{i}-1)")

summation_str = " + ".join(summation_str_parts)

# If there's more than one term in the summation, enclose it in parentheses.
if len(summation_str_parts) > 1:
    summation_str = f"({summation_str})"

# If the summation is empty (when k=n), represent it as '0'.
if not summation_str:
    summation_str = "0"

# Construct the final equation string with all numbers.
final_equation_str = f"{k}*100 - (2*{k}-1) * ({m} - {summation_str}) = {water_left}"

# Print the final resulting equation.
print(final_equation_str)

# Final answer in the required format
print(f"<<<{water_left}>>>")