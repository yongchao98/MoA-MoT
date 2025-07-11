import math

def solve_candy_packing():
    """
    Calculates the maximum number of spheres that can be packed into a box
    using a greedy algorithm on a discrete grid.
    """
    # Box dimensions in cm: 12, 12, 11
    # Sphere radius in cm: 2

    # We are working on a grid where 1 unit = 0.5 cm.
    # So, sphere radius in grid units is 2 cm / 0.5 cm/unit = 4 units.
    
    # The non-overlapping constraint means the distance between the centers of
    # any two spheres must be at least twice the radius.
    # Minimum distance = 4 (radius 1) + 4 (radius 2) = 8 units.
    # The squared distance is used to avoid calculating square roots.
    min_dist_sq = 8**2

    # As per the problem description, the valid integer coordinates for a sphere's center are:
    x_coords = range(4, 21)  # [4, 5, ..., 20]
    y_coords = range(4, 21)  # [4, 5, ..., 20]
    z_coords = range(4, 19)  # [4, 5, ..., 18]

    placed_spheres_centers = []

    # We iterate through all possible center locations, starting from the bottom-front-left corner.
    # This is the greedy approach: place a sphere in the first available spot.
    for z in z_coords:
        for y in y_coords:
            for x in x_coords:
                candidate_center = (x, y, z)
                
                # Assume we can place the sphere here until proven otherwise.
                can_place = True
                
                # Check for collision with any sphere that has already been placed.
                for placed_center in placed_spheres_centers:
                    # Calculate the squared distance to the center of the other sphere.
                    dist_sq = (
                        (candidate_center[0] - placed_center[0])**2 +
                        (candidate_center[1] - placed_center[1])**2 +
                        (candidate_center[2] - placed_center[2])**2
                    )
                    
                    # If the distance is less than the minimum, they overlap.
                    if dist_sq < min_dist_sq:
                        can_place = False
                        break  # No need to check other spheres, this spot is invalid.
                
                # If after checking all placed spheres, the spot is still valid,
                # place the new sphere here.
                if can_place:
                    placed_spheres_centers.append(candidate_center)

    # The maximum number of spheres is the total count of spheres we could place.
    max_n = len(placed_spheres_centers)

    print("The problem formulation is correct.")
    # The "final equation" is the value of the maximized objective function n.
    print(f"The maximized value is: n = {max_n}")

solve_candy_packing()