import math

def find_sphere_distances():
    """
    Calculates and prints the set of possible normalized distances (r) between
    hard spheres on a 2D plane for r <= 3.

    The distances are derived from the geometric arrangements possible in a jammed
    amorphous packing, which includes local motifs of both triangular and square
    lattices.
    """
    # These values represent the unique squared distances (r^2) derived from
    # combining possible geometries (triangular and square packings) for r <= 3.
    # r^2 = n^2+m^2 (square) or r^2 = n^2+nm+m^2 (triangular)
    squared_distances = [1, 2, 3, 4, 5, 7, 8, 9]

    print("The set of possible normalized distances r for r <= 3 are:")

    for val in squared_distances:
        # Calculate the distance r
        distance = math.sqrt(val)
        # Print the equation showing the calculation and the final value
        # formatted to two decimal places.
        print(f"Distance r = sqrt({val}) = {distance:.2f}")

if __name__ == "__main__":
    find_sphere_distances()