import math

def main():
    """
    Calculates the optimal scanner placement and coverage ratio for the Isis pyramid.
    """

    # --- Pyramid and Scanner Parameters ---
    BASE_SIDE = 150.0  # meters
    HEIGHT = 110.0     # meters
    R_LONG = 20.0      # meters
    R_SHORT = 7.0      # meters

    # --- Helper Classes and Functions ---
    class Sphere:
        def __init__(self, center, radius):
            self.center = tuple(float(c) for c in center)
            self.r = float(radius)
            self.xc, self.yc, self.zc = self.center

    def get_pyramid_half_width(z):
        """Calculates the pyramid's half-width at a given height z."""
        if 0 <= z <= HEIGHT:
            return (BASE_SIDE / 2.0) * (1.0 - z / HEIGHT)
        return 0

    def is_sphere_valid(sphere, existing_spheres):
        """Checks if a new sphere is valid (inside pyramid and not overlapping)."""
        # 1. Check if sphere is within pyramid boundaries
        # A simple, strict check: center must be in a smaller, inner pyramid.
        if not (sphere.r <= sphere.zc <= HEIGHT):
             return False
        if max(abs(sphere.xc), abs(sphere.yc)) + sphere.r > get_pyramid_half_width(sphere.zc):
            return False

        # 2. Check for overlap with existing spheres
        for s in existing_spheres:
            dist_sq = (sphere.xc - s.xc)**2 + (sphere.yc - s.yc)**2 + (sphere.zc - s.zc)**2
            min_dist = sphere.r + s.r
            # Use a small tolerance for floating point comparisons
            if dist_sq < min_dist**2 - 1e-9:
                return False # Overlap detected
        return True

    # --- Define Scanner Locations based on the plan ---
    long_range_centers = [
        # Layer 1 at z=20m (3x3 grid)
        (0, 0, 20), (40, 0, 20), (-40, 0, 20), (0, 40, 20), (0, -40, 20),
        (40, 40, 20), (40, -40, 20), (-40, 40, 20), (-40, -40, 20),
        # Layer 2 at z=60m (1 sphere)
        (0, 0, 60)
    ]

    short_range_centers = [
        # In hollows of z=20m layer
        (20, 20, 20), (20, -20, 20), (-20, 20, 20), (-20, -20, 20),
        # On top of the z=60m sphere
        (0, 0, 87),
        # Ring around z=60m sphere
        (14, 0, 60), (-14, 0, 60), (0, 14, 60), (0, -14, 60),
        (14, 14, 60), (14, -14, 60), (-14, 14, 60), (-14, -14, 60)
    ]

    # --- Verify locations and count ---
    placed_spheres = []
    n_long = 0
    m_short = 0

    for center in long_range_centers:
        s = Sphere(center, R_LONG)
        if is_sphere_valid(s, placed_spheres):
            placed_spheres.append(s)
            n_long += 1

    for center in short_range_centers:
        s = Sphere(center, R_SHORT)
        if is_sphere_valid(s, placed_spheres):
            placed_spheres.append(s)
            m_short += 1

    # --- Calculate Volumes and Coverage Ratio ---
    pyramid_volume = (1.0/3.0) * BASE_SIDE**2 * HEIGHT
    volume_long_single = (4.0/3.0) * math.pi * R_LONG**3
    volume_short_single = (4.0/3.0) * math.pi * R_SHORT**3

    total_scanned_volume = (n_long * volume_long_single) + (m_short * volume_short_single)
    coverage_ratio = (total_scanned_volume / pyramid_volume) * 100.0

    # --- Print Final Answer ---
    print("Optimal scanning configuration determined:")
    print(f"Number of long-range scans (n): {n_long}")
    print(f"Number of short-range scans (m): {m_short}")
    print("\nCalculating coverage ratio (p):")
    print(f"p = (n * V_long + m * V_short) / V_pyramid")
    print(f"p = ({n_long} * (4/3)*pi*{R_LONG:.0f}^3 + {m_short} * (4/3)*pi*{R_SHORT:.0f}^3) / ((1/3)*{BASE_SIDE:.0f}^2*{HEIGHT:.0f})")
    print(f"p = ({n_long * volume_long_single:.1f} + {m_short * volume_short_single:.1f}) / {pyramid_volume:.1f}")
    print(f"p = {total_scanned_volume:.1f} / {pyramid_volume:.1f}")
    print(f"p = {coverage_ratio:.3f}%")
    print(f"\nFinal Answer (n:m:p rounded to 0.1%):")
    print(f"{n_long}:{m_short}:{coverage_ratio:.1f}")

if __name__ == "__main__":
    main()
<<<10:13:42.9>>>