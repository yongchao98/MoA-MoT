import numpy as np

def to_plucker(point_on_line, direction_vector):
    """
    Computes the Plücker coordinates (direction, moment) for a 3D line.
    The line is defined by a point and a direction vector.
    """
    # Ensure the direction vector is a unit vector
    d = direction_vector / np.linalg.norm(direction_vector)
    # The moment vector is the cross product of the point and the direction
    m = np.cross(point_on_line, d)
    return d, m

def reciprocal_product(d1, m1, d2, m2):
    """
    Calculates the reciprocal product of two Plücker lines L1=(d1, m1) and L2=(d2, m2).
    The lines intersect if and only if this product is zero.
    """
    return np.dot(d1, m2) + np.dot(d2, m1)

# --- Simulation Setup ---
# All coordinates are in a common world reference frame.

# 1. Define a 3D point and two camera centers.
P_world = np.array([2.0, 3.0, 10.0])
C1 = np.array([0.0, 0.0, 0.0])  # Camera 1 is at the world origin
C2 = np.array([5.0, 0.0, 0.0])  # Camera 2 is on the x-axis

print("Triangulation using Plücker Coordinates Demonstration")
print("="*50)
print(f"Goal: Triangulate 3D Point P = {P_world}")
print(f"From Camera 1 (C1) at {C1} and Camera 2 (C2) at {C2}\n")

# --- Case 1: Ideal, Noise-Free Measurements ---
print("--- 1. Ideal Case (Perfect Measurements) ---")
# The direction of the ray is from the camera center to the 3D point
dir1_ideal = P_world - C1
dir2_ideal = P_world - C2

# Get Plücker coordinates for the two ideal lines
d1_ideal, m1_ideal = to_plucker(C1, dir1_ideal)
d2_ideal, m2_ideal = to_plucker(C2, dir2_ideal)

# The condition for intersection is: d1.m2 + d2.m1 = 0
side_product_ideal = reciprocal_product(d1_ideal, m1_ideal, d2_ideal, m2_ideal)

print("The two ideal rays should intersect at point P.")
print("We check this using the reciprocal product, which must be 0 for intersection.")
print(f"Reciprocal Product = d1.m2 + d2.m1")
print(f"d1 = {d1_ideal}")
print(f"m1 = {m1_ideal}")
print(f"d2 = {d2_ideal}")
print(f"m2 = {m2_ideal}")
print(f"Result = {np.dot(d1_ideal, m2_ideal):.6f} + {np.dot(d2_ideal, m1_ideal):.6f} = {side_product_ideal:.6f}")
print("Conclusion: The result is effectively zero, so the ideal lines intersect.\n")


# --- Case 2: Realistic, Noisy Measurements ---
print("--- 2. Realistic Case (With Measurement Noise) ---")
# Simulate a small error in the measurement from Camera 2.
# This results in a slightly perturbed direction vector for the second ray.
noise = np.array([0.0, 0.02, -0.01]) # A small, realistic error
dir2_noisy = (P_world + noise) - C2

# Get Plücker coordinates for the noisy line
d2_noisy, m2_noisy = to_plucker(C2, dir2_noisy)

# Check the intersection condition again with the noisy line
side_product_noisy = reciprocal_product(d1_ideal, m1_ideal, d2_noisy, m2_noisy)

print("Due to measurement noise, the second ray's direction is slightly off.")
print("We check the reciprocal product again.")
print(f"Reciprocal Product = d1.m2_noisy + d2_noisy.m1")
print(f"d1 = {d1_ideal}")
print(f"m1 = {m1_ideal}")
print(f"d2_noisy = {d2_noisy}")
print(f"m2_noisy = {m2_noisy}")
print(f"Result = {np.dot(d1_ideal, m2_noisy):.6f} + {np.dot(d2_noisy, m1_ideal):.6f} = {side_product_noisy:.6f}")
print("\nConclusion: The result is non-zero. The lines are skew and DO NOT intersect.")
print("\nBecause the lines do not intersect in a realistic scenario, a direct solution for an intersection point is not possible. This is a fundamental limitation.")
