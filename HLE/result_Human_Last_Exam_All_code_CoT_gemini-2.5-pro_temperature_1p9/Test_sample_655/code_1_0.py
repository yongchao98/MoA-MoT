import math

def solve_pioneer_probe_packing():
    """
    This script finds the optimal container and the number of energy balls 
    to maximize the total energy for the Pioneer probe mission.
    It works by running a greedy packing simulation on the most promising
    container shapes that respect the surface area constraint.
    """

    class Container:
        """Base class for container shapes."""
        def get_bounds(self):
            raise NotImplementedError
        def is_inside(self, point, radius):
            raise NotImplementedError
        def description(self):
            raise NotImplementedError
        @property
        def surface_area(self):
            raise NotImplementedError
        @property
        def volume(self):
            raise NotImplementedError

    class Box(Container):
        """Represents a box-shaped container aligned with axes from origin."""
        def __init__(self, l, w, h):
            self.l, self.w, self.h = l, w, h

        def get_bounds(self):
            return (0, self.l, 0, self.w, 0, self.h)

        def is_inside(self, point, radius):
            x, y, z = point
            # Use a small tolerance for floating point comparisons
            tol = 1e-9
            return (radius - tol <= x <= self.l - radius + tol and
                    radius - tol <= y <= self.w - radius + tol and
                    radius - tol <= z <= self.h - radius + tol)

        def description(self):
            return f"box {self.l}x{self.w}x{self.h}"

        @property
        def surface_area(self):
            return 2 * (self.l * self.w + self.w * self.h + self.h * self.l)
        
        @property
        def volume(self):
            return self.l * self.w * self.h

    class Cylinder(Container):
        """Represents a cylindrical container with its base at z=0 and center on the z-axis."""
        def __init__(self, r, h):
            self.r, self.h = r, h

        def get_bounds(self):
            return (-self.r, self.r, -self.r, self.r, 0, self.h)

        def is_inside(self, point, radius):
            x, y, z = point
            tol = 1e-9
            return (math.sqrt(x**2 + y**2) + radius <= self.r + tol and
                    radius - tol <= z <= self.h - radius + tol)

        def description(self):
            return f"cylinder r={self.r}, h={self.h}"

        @property
        def surface_area(self):
            return 2 * math.pi * self.r**2 + 2 * math.pi * self.r * self.h
        
        @property
        def volume(self):
            return math.pi * self.r**2 * self.h

    class Sphere(Container):
        """Represents a spherical container centered at the origin."""
        def __init__(self, r):
            self.r = r

        def get_bounds(self):
            return (-self.r, self.r, -self.r, self.r, -self.r, self.r)

        def is_inside(self, point, radius):
            tol = 1e-9
            return math.sqrt(point[0]**2 + point[1]**2 + point[2]**2) + radius <= self.r + tol

        def description(self):
            return f"sphere r={self.r}"

        @property
        def surface_area(self):
            return 4 * math.pi * self.r**2
        
        @property
        def volume(self):
            return 4/3 * math.pi * self.r**3

    def pack_container(container, grid_step=0.5):
        """Greedy algorithm to pack spheres into a container."""
        placed_balls = []
        
        # 1. Pack large balls (radius=2)
        r_large = 2.0
        min_x, max_x, min_y, max_y, min_z, max_z = container.get_bounds()
        candidate_centers_large = []
        for kx in range(math.ceil(min_x / grid_step), math.floor(max_x / grid_step) + 1):
            for ky in range(math.ceil(min_y / grid_step), math.floor(max_y / grid_step) + 1):
                for kz in range(math.ceil(min_z / grid_step), math.floor(max_z / grid_step) + 1):
                    p = (kx * grid_step, ky * grid_step, kz * grid_step)
                    if container.is_inside(p, r_large):
                        candidate_centers_large.append(p)
        
        candidate_centers_large.sort(key=lambda p: (p[2], p[1], p[0]))

        for p in candidate_centers_large:
            can_place = True
            for ball in placed_balls:
                dist_sq = sum((p[i] - ball['center'][i])**2 for i in range(3))
                if dist_sq < (r_large + ball['radius'])**2 - 1e-9:
                    can_place = False
                    break
            if can_place:
                placed_balls.append({'center': p, 'radius': r_large})
        n2 = len(placed_balls)

        # 2. Pack small balls (radius=1)
        r_small = 1.0
        candidate_centers_small = []
        for kx in range(math.ceil(min_x / grid_step), math.floor(max_x / grid_step) + 1):
            for ky in range(math.ceil(min_y / grid_step), math.floor(max_y / grid_step) + 1):
                for kz in range(math.ceil(min_z / grid_step), math.floor(max_z / grid_step) + 1):
                    p = (kx * grid_step, ky * grid_step, kz * grid_step)
                    if container.is_inside(p, r_small):
                         is_large_ball_center = any(all(abs(p[i] - ball['center'][i]) < 1e-9 for i in range(3)) for ball in placed_balls if ball['radius'] == r_large)
                         if not is_large_ball_center:
                             candidate_centers_small.append(p)
        
        candidate_centers_small.sort(key=lambda p: (p[2], p[1], p[0]))

        for p in candidate_centers_small:
            can_place = True
            for ball in placed_balls:
                dist_sq = sum((p[i] - ball['center'][i])**2 for i in range(3))
                if dist_sq < (r_small + ball['radius'])**2 - 1e-9:
                    can_place = False
                    break
            if can_place:
                placed_balls.append({'center': p, 'radius': r_small})
                
        n1 = len(placed_balls) - n2
        energy = n1 * 1 + n2 * 10
        return energy, n1, n2

    # --- Main Execution ---
    print("Starting simulation to find the best container for Pioneer probe...")
    # These candidates are optimized for maximum volume under the surface area constraint of 1050 cm^2
    candidates = [
        Sphere(9.0),
        Cylinder(7.5, 14.5),
        Box(13.0, 13.0, 13.0),
    ]

    max_energy = -1
    best_result = None
    max_surface_area = 1050.0

    for container in candidates:
        if container.surface_area > max_surface_area + 1e-9:
            continue
        
        print(f"\nTesting {container.description()} (SA={container.surface_area:.1f} cm^2, Vol={container.volume:.1f} cm^3)")
        energy, n1, n2 = pack_container(container)
        print(f" -> Result: {n2} large balls, {n1} small balls. Total Energy: {energy} MJ.")

        if energy > max_energy:
            max_energy = energy
            best_result = (container, n1, n2)

    print("\n--- Optimal Configuration Found ---")
    if best_result:
        container, n1, n2 = best_result
        final_answer_str = f"[{container.description()}]a={n1};b={n2}"
        print(final_answer_str)
        print("Final equation for energy calculation:")
        print(f"1 * {n1} + 10 * {n2} = {max_energy}")
        print(f"\n<<<{final_answer_str}>>>")
    else:
        print("No valid solution was found.")

solve_pioneer_probe_packing()