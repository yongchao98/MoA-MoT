import math

def calculate_areas():
    """
    Calculates the area of various shapes with a fixed perimeter of 1 meter.
    This demonstrates the isoperimetric inequality principle.
    """
    perimeter = 1.0
    
    # A. Equilateral Triangle
    side_triangle = perimeter / 3
    area_triangle = (math.sqrt(3) / 4) * side_triangle**2
    
    # D. Square
    side_square = perimeter / 4
    area_square = side_square**2
    
    # E. Regular Hexagon
    side_hexagon = perimeter / 6
    area_hexagon = (3 * math.sqrt(3) / 2) * side_hexagon**2
    
    # G. Circle
    radius_circle = perimeter / (2 * math.pi)
    area_circle = math.pi * radius_circle**2
    
    # Store results for comparison
    results = {
        "Equilateral Triangle": area_triangle,
        "Square": area_square,
        "Regular Hexagon": area_hexagon,
        "Circle": area_circle,
    }
    
    print(f"For a fixed perimeter of {perimeter} meter:\n")
    for shape, area in results.items():
        print(f"Area of {shape}: {area:.6f} m^2")
        
    # Find the shape with the maximum area
    max_area_shape = max(results, key=results.get)
    
    print(f"\nConclusion: The {max_area_shape} encloses the maximum area.")
    print("Therefore, the cut structure which optimizes the area has the symmetry group of a circle.")

if __name__ == '__main__':
    calculate_areas()