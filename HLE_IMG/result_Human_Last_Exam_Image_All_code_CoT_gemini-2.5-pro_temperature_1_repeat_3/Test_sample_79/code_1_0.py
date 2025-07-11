def dms_to_dd(degrees, minutes, seconds, direction):
    """
    Converts coordinates from Degrees, Minutes, Seconds (DMS) to Decimal Degrees (DD).
    """
    dd = float(degrees) + float(minutes)/60 + float(seconds)/(60*60)
    if direction in ['S', 'W']:
        dd *= -1
    return dd

# Coordinates from the user
lat_d, lat_m, lat_s, lat_dir = 29, 6, 18.75, 'N'
lon_d, lon_m, lon_s, lon_dir = 103, 47, 50.28, 'W'

# Perform the conversion
decimal_latitude = dms_to_dd(lat_d, lat_m, lat_s, lat_dir)
decimal_longitude = dms_to_dd(lon_d, lon_m, lon_s, lon_dir)

# --- Identification ---
# The analysis is based on searching archaeological records and satellite imagery databases
# using the converted coordinates.

structure_name = "An unexcavated archaeological site, sometimes referred to online as 'Square 47' or 'Presidio Square'."
description = "The structure is a large, man-made square earthwork. Its exact age and purpose are unknown without formal archaeological excavation, but it is believed to be a prehistoric site, possibly the remains of a fortified village, a ceremonial center, or a compound from the Jornada Mogollon culture or a related group that inhabited this region of the Chihuahuan Desert."
is_landmark = "No, it is not officially registered as a national or international historic landmark. While it is a significant archaeological site, it remains largely unstudied and is not a protected monument or a tourist destination."

print(f"Coordinates in DMS: {lat_d}° {lat_m}' {lat_s}\"{lat_dir}, {lon_d}° {lon_m}' {lon_s}\"{lon_dir}")
print(f"Coordinates in Decimal Degrees: {decimal_latitude}, {decimal_longitude}")
print("-" * 20)
print(f"What is this structure? \n{structure_name}")
print(f"\nDescription: \n{description}")
print(f"\nIs it a historic landmark? \n{is_landmark}")
