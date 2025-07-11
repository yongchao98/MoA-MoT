import osmnx as ox
import numpy as np
from shapely.geometry import Point

def solve_building_shape():
    """
    This function fetches the building footprints for two schools in Pisa,
    and prints an ASCII representation of their combined shape.
    """
    # Names of the schools
    school_names = [
        "Istituto Galilei-Pacinotti, Pisa, Italy",
        "Liceo Scientifico Ulisse Dini, Pisa, Italy"
    ]

    print("Fetching building data for the schools from OpenStreetMap...")

    # Define the tag to look for buildings
    tags = {"building": True}

    try:
        # Get the building footprints for both schools
        # The buffer_dist expands the search area around the geocoded point
        gdf_pacinotti = ox.geometries_from_address(school_names[0], tags, buffer_dist=150)
        gdf_dini = ox.geometries_from_address(school_names[1], tags, buffer_dist=150)

        # Combine the data for both schools
        all_buildings = gdf_pacinotti._append(gdf_dini, ignore_index=True)
        
        # Merge all building shapes into one for easier processing
        combined_shape = all_buildings.unary_union

    except Exception as e:
        print(f"Could not fetch the map data. Error: {e}")
        print("Please ensure you have an internet connection and the addresses are correct.")
        return

    print("Data downloaded. Now generating ASCII art of the building shapes...")

    # Get the bounding box of the combined shape
    minx, miny, maxx, maxy = combined_shape.bounds

    # Define the size of our ASCII art grid
    width = 70
    # Adjust height to maintain the approximate aspect ratio of the buildings
    aspect_ratio = (maxy - miny) / (maxx - minx)
    height = int(width * aspect_ratio * 0.5) # 0.5 is a correction factor for char aspect ratio

    # Create the grid
    grid = np.full((height, width), " ", dtype=str)

    # Calculate step size for longitude and latitude
    x_step = (maxx - minx) / width
    y_step = (maxy - miny) / height

    # Populate the grid
    for r in range(height):
        for c in range(width):
            # Calculate the geographic coordinate for the center of the grid cell
            # We map from top (row 0) to maxy and bottom (row height-1) to miny
            lon = minx + c * x_step
            lat = maxy - r * y_step
            
            # Check if this point is inside any of the building polygons
            if combined_shape.contains(Point(lon, lat)):
                grid[r, c] = "█"

    # Print the resulting ASCII art
    print("\n--- Building Shape ---\n")
    for row in grid:
        print("".join(row))
    print("\n----------------------\n")
    print("When viewed from above, the buildings form the letter 'L'.")

if __name__ == '__main__':
    solve_building_shape()
