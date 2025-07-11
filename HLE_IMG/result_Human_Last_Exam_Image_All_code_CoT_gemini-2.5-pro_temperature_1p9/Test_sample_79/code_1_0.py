import sys

def identify_structure():
    """
    Identifies the structure at the given coordinates and provides information.
    The coordinates point to a well-known land art installation.
    """
    # Coordinates provided by the user
    lat_deg = 29
    lat_min = 6
    lat_sec = 18.75
    # The given latitude is North (positive)

    lon_deg = 103
    lon_min = 47
    lon_sec = 50.28
    # The given longitude is West (negative)

    # --- Coordinate Conversion Calculation ---
    # The prompt requires showing the numbers in the equation.
    print("--- Coordinate Conversion ---\n")
    
    # Calculate latitude in decimal degrees
    decimal_lat = lat_deg + (lat_min / 60) + (lat_sec / 3600)
    print("Calculating Decimal Latitude (N):")
    print(f"{lat_deg} + ({lat_min} / 60) + ({lat_sec} / 3600) = {decimal_lat:.6f}\n")

    # Calculate longitude in decimal degrees
    decimal_lon = -1 * (lon_deg + (lon_min / 60) + (lon_sec / 3600))
    print("Calculating Decimal Longitude (W):")
    print(f"-1 * ({lon_deg} + ({lon_min} / 60) + ({lon_sec} / 3600)) = {decimal_lon:.6f}\n")

    # --- Structure Identification ---
    print("--- Structure Information ---\n")
    print("The feature at these coordinates is NOT an ancient structure or historic ruin.\n")
    print("Identification: It is '15 untitled works in concrete' by the artist Donald Judd.")
    print("Description: This is a large-scale contemporary land art installation created between 1980-1984.")
    print("The work consists of fifteen large, hollow concrete boxes arranged in a straight line extending one kilometer through the Chihuahuan Desert.\n")
    
    # --- Landmark Status ---
    print("--- Landmark Status ---\n")
    print("While not a traditional 'historic' landmark, it is a globally significant cultural landmark of modern and minimalist art.")
    print("It is managed by the Chinati Foundation, a major art institution based in Marfa, Texas, and is a key destination for art enthusiasts.")

if __name__ == "__main__":
    identify_structure()