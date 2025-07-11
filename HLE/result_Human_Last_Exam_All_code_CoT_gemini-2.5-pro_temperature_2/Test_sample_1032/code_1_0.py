import math

# Step 1: Define parameters
L = 7 / 2  # Rope length
house_area = 3.0  # Area of the house (3 unit squares)

# Step 2: Calculate the area in Quadrants 1, 2, and 4 (unobstructed)
# The area of a quarter taxi-cab circle is L^2 / 2.
area_q1_q2_q4 = 3 * (L**2) / 2

# Step 3: Calculate the directly accessible area in Quadrant 3
# This is the area of the Q3 quarter circle minus the house area.
area_q3_unconstrained = (L**2) / 2
area_q3_direct = area_q3_unconstrained - house_area

# Step 4: Calculate the "shadow" areas accessible by the rope bending around corners.
# The calculation for the total shadow area on each side (e.g., left, bottom)
# simplifies to the area of the quarter-circle formed from the outermost pivot point.

# Shadow region to the left of the house (x < -2).
# The pivot is at C_x = (-2, 0). The rope distance to this pivot is 2.
remaining_rope_left = L - 2
shadow_area_left = (remaining_rope_left**2) / 2

# Shadow region below the house (y < -2).
# The pivot is at C_y = (0, -2). The rope distance to this pivot is 2.
remaining_rope_bottom = L - 2
shadow_area_bottom = (remaining_rope_bottom**2) / 2

# Shadow region in the middle pocket (-2 < x < -1, -2 < y < -1).
# This area is the union of two disjoint quarter-circles from two inner pivots.
# Pivot 1 is C1(-2, -1), reached via path of length 3. Remaining rope L - 3 = 0.5.
remaining_rope_C1 = L - 3
shadow_area_C1 = (remaining_rope_C1**2) / 2
# Pivot 2 is C2(-1, -2), reached via path of length 3. Remaining rope L - 3 = 0.5.
remaining_rope_C2 = L - 3
shadow_area_C2 = (remaining_rope_C2**2) / 2
shadow_area_middle = shadow_area_C1 + shadow_area_C2

# Step 5: Sum all the disjoint component areas for the final result.
total_area = area_q1_q2_q4 + area_q3_direct + shadow_area_left + shadow_area_bottom + shadow_area_middle

# Step 6: Print the breakdown of the final calculation as requested.
print("The total area is the sum of the following disjoint regions:")
print(f"1. Unobstructed Area in Q1, Q2, Q4: {area_q1_q2_q4}")
print(f"2. Directly Accessible Area in Q3: {area_q3_direct}")
print(f"3. Shadow Area Left of House (x<-2): {shadow_area_left}")
print(f"4. Shadow Area Below House (y<-2): {shadow_area_bottom}")
print(f"5. Shadow Area in Middle Pocket: {shadow_area_middle}")
print("\nFinal Equation:")
print(f"Total Area = {area_q1_q2_q4} + {area_q3_direct} + {shadow_area_left} + {shadow_area_bottom} + {shadow_area_middle}")
print(f"Total Area = {total_area}")
