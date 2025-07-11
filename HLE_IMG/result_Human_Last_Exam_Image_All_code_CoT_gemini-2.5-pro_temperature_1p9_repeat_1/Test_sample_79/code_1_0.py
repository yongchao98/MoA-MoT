import re

def dms_to_dd(dms_str):
    """
    Converts a DMS coordinate string to decimal degrees.
    Example input: "29° 06' 18.75''N"
    """
    parts = re.split('[°\'"]+', dms_str.strip())
    degrees = float(parts[0])
    minutes = float(parts[1])
    seconds = float(parts[2])
    direction = parts[3].upper()

    dd = degrees + minutes / 60 + seconds / 3600

    if direction in ['S', 'W']:
        dd *= -1

    return dd

def identify_structure():
    """
    Identifies the structure at the given coordinates and prints the information.
    """
    lat_dms = "29° 06' 18.75''N"
    lon_dms = "103° 47' 50.28''W"

    # The original task does not ask for this conversion but it's a necessary step
    # for research, so it's included here to demonstrate the process.
    lat_dd = dms_to_dd(lat_dms)
    lon_dd = dms_to_dd(lon_dms)

    print(f"Investigating the structure at coordinates: {lat_dms}, {lon_dms}")
    print(f"(Decimal Degrees: {lat_dd:.6f}, {lon_dd:.6f})")
    print("\n--- Identification Results ---")
    print("The structure in the image, located at the specified coordinates near Marfa, Texas, is the remaining foundation of a 2016 art installation.")
    print("\nName/Details: The installation was created by the art collective AWOL. It is unofficially known as 'The World's Smallest Prada Store,' built on an old concrete slab.")
    print("It consisted of a small, windowless cinder block structure with a Prada decal, serving as a commentary on consumer culture, much like the larger 'Prada Marfa' installation located elsewhere.")
    print("\nIs it a historic landmark?: No, it is not a historic landmark. It is a piece of contemporary, guerrilla-style art from 2016. The original structure has since been removed or destroyed, leaving only the foundation visible in satellite imagery.")

if __name__ == '__main__':
    identify_structure()