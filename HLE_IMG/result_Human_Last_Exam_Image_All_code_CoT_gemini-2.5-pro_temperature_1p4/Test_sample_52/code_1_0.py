import math
from fractions import Fraction

def main():
    """
    Calculates the tiling period for four given shapes on a hexagonal grid.
    """
    # Step 1: Define the coordinate system and point locations.
    # Hexagons are flat-top with side length a=2. s3 is sqrt(3).
    s3 = math.sqrt(3)
    
    # Coordinates of all 31 points
    coords = {
        # H_13 (center 13) perimeter and center
        1: (-2, 0), 2: (-1.5, -s3/2), 3: (-1, -s3), 4: (0, -s3), 5: (1, -s3),
        6: (1.5, -s3/2), 7: (2, 0), 8: (1.5, s3/2), 9: (1, s3), 10: (0, s3),
        11: (-1, s3), 12: (-1.5, s3/2), 13: (0, 0),
        # H_23 (center 23) perimeter and center
        14: (1.5, -1.5*s3), 15: (2, -2*s3), 16: (3, -2*s3), 17: (4, -2*s3),
        18: (4.5, -1.5*s3), 19: (5, -s3), 20: (4.5, -s3/2), 21: (4, 0),
        22: (3, 0), 23: (3, -s3),
        # H_31 (center 31) perimeter and center
        24: (4.5, s3/2), 25: (5, s3), 26: (4.5, 1.5*s3), 27: (4, 2*s3),
        28: (3, 2*s3), 29: (2, 2*s3), 30: (1.5, 1.5*s3), 31: (3, s3)
    }

    # Step 2: Define helper functions for area and period calculation.
    def polygon_area(points):
        """Calculates polygon area using the shoelace formula."""
        area = 0.0
        for i in range(len(points)):
            p1 = points[i]
            p2 = points[(i + 1) % len(points)]
            area += p1[0] * p2[1] - p2[0] * p1[1]
        return abs(area) / 2.0

    # Step 3: Process each of the four cases.
    cases = [
        [13, 31, 23],
        [10, 4, 23, 31],
        [5, 15, 17, 19, 21, 7],
        [4, 5, 14, 23, 18, 19, 21, 22, 31, 30, 9, 8, 13]
    ]

    periods = []
    a_hex_scaled = 6  # Area of one hexagon is 6 * sqrt(3)

    for i, case_points_labels in enumerate(cases):
        print(f"--- Case {i+1} ---")
        points_str = ", ".join(map(str, case_points_labels))
        print(f"Points defining the tile: {points_str}")

        case_coords = [coords[p] for p in case_points_labels]
        
        # Calculate tile area and express it as a rational multiple of sqrt(3)
        tile_area = polygon_area(case_coords)
        tile_area_scaled = tile_area / s3
        
        frac = Fraction(tile_area_scaled).limit_denominator(1000)
        area_num = frac.numerator
        area_den = frac.denominator
        
        print(f"Tile Area = {tile_area:.4f}")
        print(f"This area is equivalent to ({area_num}/{area_den}) * sqrt(3).")
        
        # Calculate period 'k' using the area relation: k * A_tile = n * A_hex
        # k * (area_num/area_den) * s3 = n * a_hex_scaled * s3
        # k * area_num = n * a_hex_scaled * area_den
        # We need the smallest integer k, which means we simplify the ratio.
        g = math.gcd(area_num, a_hex_scaled * area_den)
        k = (a_hex_scaled * area_den) // g
        n = area_num // g
        
        print(f"The tiling condition is given by the equation: k * (Area of Tile) = n * (Area of Hexagon).")
        print(f"Substituting the area values (scaled by sqrt(3)):")
        print(f"k * {area_num}/{area_den} = n * {a_hex_scaled}")
        print(f"This simplifies to the Diophantine equation: k * {area_num} = n * {a_hex_scaled * area_den}")
        print(f"The smallest integer solution is k = {k} and n = {n}.")
        print(f"The period is {k}.\n")
        periods.append(k)

    # Step 4: Output the final combined answer.
    final_answer = ",".join(map(str, periods))
    print(f"--- Final Answer ---")
    print(f"The four periods separated by commas are: {final_answer}")
    print(f"<<<{final_answer}>>>")

if __name__ == '__main__':
    main()