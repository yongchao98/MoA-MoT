def dms_to_decimal(degrees, minutes, seconds):
    """Converts coordinates from DMS to decimal degrees."""
    return degrees + minutes / 60 + seconds / 3600

# --- Latitude Calculation ---
lat_deg = 29
lat_min = 6
lat_sec = 18.75

# --- Longitude Calculation ---
lon_deg = 103
lon_min = 47
lon_sec = 50.28

# Perform the conversions
decimal_lat = dms_to_decimal(lat_deg, lat_min, lat_sec)
# Longitude is West, so it's negative
decimal_lon = -dms_to_decimal(lon_deg, lon_min, lon_sec)

print("Step 1: Convert coordinates to decimal format for identification.\n")

print("Latitude Conversion:")
print(f"{lat_deg}° + {lat_min}'/60 + {lat_sec}\"/3600 = {decimal_lat:.6f}° N")

print("\nLongitude Conversion:")
# Make the final result negative for the printout
print(f"-({lon_deg}° + {lon_min}'/60 + {lon_sec}\"/3600) = {decimal_lon:.6f}° W")

print("\n--------------------------------------------------\n")
print("Step 2: Identify the structure at the location.\n")

print("What is the structure?")
print("The square object at these coordinates is a military aerial photography calibration target.")
print("It is known as the 'Marfa Army Airfield Calibration Target' or 'Marfa Grid'.")
print("These large, concrete targets were built by the U.S. Army Air Forces during and after World War II.")
print("Their purpose was to allow high-altitude reconnaissance aircraft to test and calibrate the resolution and focus of their cameras.\n")

print("Is it a historic landmark?")
print("While it is a significant historical relic from the Cold War era, it is not an officially designated national or state historic landmark. It is a point of interest for aviation and military history enthusiasts.")
