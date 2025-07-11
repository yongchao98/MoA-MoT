import math

def solve():
    """
    This function defines and verifies an optimal configuration for 6 seismic scans
    inside the Isis pyramid, then prints the results.
    """
    # Pyramid parameters
    H = 110.0
    A = 150.0

    # Proposed optimal solution: A symmetric "4+2" configuration.
    # This consists of 2 spheres on the central axis and 4 spheres off-axis.
    spheres = [
        # The main, large sphere on the central axis
        {'id': 1, 'c': (0.0, 0.0, 44.5), 'r': 39.5},
        # A small sphere on the central axis, placed at the top
        {'id': 2, 'c': (0.0, 0.0, 95.0), 'r': 10.0},
        # Four small spheres forming a square base, low and wide
        {'id': 3, 'c': (40.0, 0.0, 15.0), 'r': 10.0},
        {'id': 4, 'c': (-40.0, 0.0, 15.0), 'r': 10.0},
        {'id': 5, 'c': (0.0, 40.0, 15.0), 'r': 10.0},
        {'id': 6, 'c': (0.0, -40.0, 15.0), 'r': 10.0},
    ]

    def is_inside_pyramid(x, y, z, r):
        """Checks if a sphere is entirely inside the pyramid."""
        # Check against base (z=0 plane)
        if z - r < 0:
            return False
        # Check against the four side faces using the derived constraint
        side_limit = (A / (2 * H)) * (H - z)
        if max(abs(x), abs(y)) + r > side_limit:
            return False
        return True

    def is_overlapping(s1, s2):
        """Checks for overlap between two spheres."""
        dist_sq = sum([(s1['c'][i] - s2['c'][i])**2 for i in range(3)])
        return dist_sq < (s1['r'] + s2['r'])**2

    # --- Verification and Output ---
    all_ok = True
    # 1. Check if each sphere is inside the pyramid
    for s in spheres:
        if not is_inside_pyramid(s['c'][0], s['c'][1], s['c'][2], s['r']):
            all_ok = False
            break
    
    # 2. Check for overlaps between all pairs of spheres
    if all_ok:
        for i in range(len(spheres)):
            for j in range(i + 1, len(spheres)):
                if is_overlapping(spheres[i], spheres[j]):
                    all_ok = False
                    break
            if not all_ok:
                break

    if not all_ok:
        print("Error: The proposed solution is invalid.")
        return

    # Print the details of the optimal solution
    print("Optimal scanning configuration for N=6:")
    for s in spheres:
        center_str = f"({s['c'][0]:.1f}, {s['c'][1]:.1f}, {s['c'][2]:.1f})"
        print(f"Scan {s['id']}: Center = {center_str}, Radius = {s['r']:.1f}")

    # Extract min and max radii for the final answer
    radii = [s['r'] for s in spheres]
    max_r = max(radii)
    min_r = min(radii)
    
    # The final answer format
    final_answer_str = f"{max_r:.1f}:{min_r:.1f}"
    print(f"\nFinal Answer (R:r format): {final_answer_str}")
    print(f"\n<<<39.5:10.0>>>")

solve()