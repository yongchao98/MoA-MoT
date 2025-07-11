import math

def dms_to_dd(degrees, minutes, seconds):
    """Converts Degrees, Minutes, Seconds to Decimal Degrees."""
    dd = float(degrees) + float(minutes)/60 + float(seconds)/3600
    return dd

# --- Input Coordinates ---
# Latitude: 29° 3' 28.15''N
lat_d = 29
lat_m = 3
lat_s = 28.15

# Longitude: 103° 48' 11.84''W
lon_d = 103
lon_m = 48
lon_s = 11.84

# --- Coordinate Conversion Calculation ---
print("To analyze the location, we convert the coordinates from Degrees/Minutes/Seconds (DMS) to Decimal Degrees (DD).")

# Perform and display latitude conversion
lat_dd = dms_to_dd(lat_d, lat_m, lat_s)
print("\n--- Latitude Conversion (29° 3' 28.15''N) ---")
print(f"Equation: Degrees + (Minutes / 60) + (Seconds / 3600)")
print(f"Calculation: {lat_d} + ({lat_m} / 60) + ({lat_s} / 3600) = {lat_dd:.6f}° N")

# Perform and display longitude conversion
# West longitude is negative in decimal degrees
lon_dd = -dms_to_dd(lon_d, lon_m, lon_s)
print("\n--- Longitude Conversion (103° 48' 11.84''W) ---")
print(f"Equation: Degrees + (Minutes / 60) + (Seconds / 3600)")
print(f"Calculation: {lon_d} + ({lon_m} / 60) + ({lon_s} / 3600) = {-lon_dd:.6f}°")
print(f"Result for mapping: {lon_dd:.6f}° W")

# --- Analysis ---
print("\n--- Analysis of Carving and Location ---")
print("1. Location: The coordinates place the rock in the Big Bend region of Texas, an area with a rich history of indigenous peoples and rock art.")
print("2. Image Analysis: The carving's most prominent feature is a circular depression made of many pecked dots. This pattern does not clearly resemble a river's path or a mountain range's silhouette.")
print("3. Evaluating Options:")
print("   - Options B and D claim the carving is a map. However, the directions and distances provided (e.g., Chisos 'north', Bravo River 'northwest') are geographically inaccurate for this location. The carving's form is also not a clear topographical match.")
print("   - Option C is false, as carvings are clearly visible.")
print("   - Option E suggests the carving is a 'sipapu'. A sipapu is a symbolic feature in Southwestern indigenous cosmology, representing the portal of emergence from the underworld. In rock art, it is often depicted as a pecked pit, circle, or spiral. The carving in the image is a strong visual match for a sipapu petroglyph.")
print("\nConclusion: The carving is not a geographical map but a symbolic or religious image.")
<<<E>>>