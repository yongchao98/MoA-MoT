import numpy as np
import math

# --- Problem Constants ---
ROOM_W = 140
ROOM_H = 110
ROOM_AREA = ROOM_W * ROOM_H
COVERAGE_TARGET_RATIO = 0.88
TARGET_COVERED_AREA = ROOM_AREA * COVERAGE_TARGET_RATIO

# --- Chosen Scanner Configuration ---
# Based on cost-effectiveness analysis, we use only the C2 scanner.
# We test a 4x3 grid arrangement, totaling 12 scanners.
NUM_C2 = 12
NUM_C1 = 0
NUM_R1 = 0

C2_COST = 20000
C2_RADIUS = 20

# Define the centers for the 12 scanners in a 4x3 grid
# These coordinates are all multiples of 5.
x_centers = [20, 55, 90, 125]
y_centers = [20, 55, 90]
scanner_centers = []
for x in x_centers:
    for y in y_centers:
        scanner_centers.append((x, y))

# --- Simulation to Verify Coverage ---
# Create a grid representing 1m x 1m plots of the room
is_covered = np.zeros((ROOM_W + 1, ROOM_H + 1), dtype=bool)

# Mark the areas covered by each of the 12 scanners
for cx, cy in scanner_centers:
    # Define the bounding box for the circle to check
    x_min, x_max = max(0, cx - C2_RADIUS), min(ROOM_W, cx + C2_RADIUS)
    y_min, y_max = max(0, cy - C2_RADIUS), min(ROOM_H, cy + C2_RADIUS)
    
    for x in range(int(x_min), int(x_max) + 1):
        for y in range(int(y_min), int(y_max) + 1):
            # If the point is within the circle's radius, mark it as covered
            if (x - cx)**2 + (y - cy)**2 <= C2_RADIUS**2:
                is_covered[x, y] = True

# Calculate total number of covered 1m x 1m plots
total_covered_points = np.sum(is_covered)
final_coverage_ratio = total_covered_points / ROOM_AREA

# --- Final Calculation and Output ---
total_cost = NUM_C2 * C2_COST + NUM_C1 * 1600 + NUM_R1 * 2000

print("--- Optimal Solution Verification ---")
print(f"Proposed configuration: {NUM_C2} C2 scanners, {NUM_C1} C1 scanners, {NUM_R1} R1 scanners.")
print(f"Simulated coverage area: {total_covered_points} m^2")
print(f"Room area: {ROOM_AREA} m^2")
print(f"Target coverage area: >{TARGET_COVERED_AREA:.0f} m^2")
print(f"Final coverage ratio: {final_coverage_ratio:.4f} (Target was {COVERAGE_TARGET_RATIO})")
if total_covered_points >= TARGET_COVERED_AREA:
    print("Coverage requirement is met.")
else:
    print("Coverage requirement is NOT met.")

print("\n--- Final Answer ---")
# The final equation is printed with each number for clarity
final_equation = f"{NUM_C2} * {C2_COST}"
print(f"Optimal total cost equation: {final_equation} = {total_cost}")
