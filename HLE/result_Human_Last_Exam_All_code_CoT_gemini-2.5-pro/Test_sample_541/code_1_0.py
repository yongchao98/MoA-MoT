import math

def calculate_hull_area(n, total_length=1.0):
    """
    Calculates the area of the convex hull of an n-pointed star.
    The convex hull is a regular n-gon.

    Args:
        n (int): The number of points/arms in the star. Must be >= 2.
        total_length (float): The total length of all arms combined.

    Returns:
        float: The area of the convex hull.
    """
    if n < 2:
        return 0.0
    
    # For a total length of 1, each of the n arms has length 1/n.
    # This length is the circumradius of the n-gon.
    r = total_length / n
    
    # The area of a regular n-gon with circumradius r is (1/2) * n * r^2 * sin(2*pi/n)
    angle = 2 * math.pi / n
    area = 0.5 * n * (r**2) * math.sin(angle)
    return area, r, angle

def main():
    """
    Compares the convex hull areas for different star-shaped cuts
    and identifies the one with the maximum area.
    """
    symmetries = {
        3: "equilateral triangle",
        4: "square",
        6: "regular hexagon"
    }
    
    max_area = 0
    best_n = 0
    
    print("Calculating the area of the convex hull for different symmetric cuts (total length = 1m):\n")

    for n in sorted(symmetries.keys()):
        area, arm_length, angle_rad = calculate_hull_area(n)
        
        if area > max_area:
            max_area = area
            best_n = n
            
        print(f"Symmetry: {symmetries[n].capitalize()} ({n}-pointed star)")
        print(f"Equation: Area = (1/2) * n * r^2 * sin(2*pi/n)")
        print(f"Calculation: Area = (1/2) * {n} * ({1.0/n:.4f})^2 * sin(2*pi/{n})")
        print(f"             = {0.5 * n:.1f} * {arm_length**2:.4f} * sin({angle_rad:.4f})")
        print(f"             = {0.5 * n:.1f} * {arm_length**2:.4f} * {math.sin(angle_rad):.4f}")
        print(f"Resulting Area: {area:.6f} m^2\n")

    print("--- Conclusion ---")
    print(f"The maximum area ({max_area:.6f} m^2) is achieved with n = {best_n}.")
    print(f"This corresponds to a 3-pointed star with the symmetry of an {symmetries[best_n]}.")

if __name__ == "__main__":
    main()