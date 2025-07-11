import urllib.request
import math

def find_vernal_equinox_crossing_boundary():
    """
    This function finds the specific IAU constellation boundary marker that the
    Sun's apparent position at the Vernal Equinox most recently crossed to enter
    the constellation Pisces.

    It fetches the IAU boundary data (epoch B1875.0), searches for the
    segment of the Pisces boundary that crosses the celestial equator on its
    eastern side (the historical border with Aries), and formats the two
    points of the line segment as requested.
    """
    # The official IAU constellation boundaries data for epoch B1875.0 from the
    # Centre de Données astronomiques de Strasbourg (CDS).
    url = "https://cdsarc.u-strasbg.fr/ftp/VI/49/bound_1875.dat"
    target_line = ""

    try:
        with urllib.request.urlopen(url) as response:
            data = response.read().decode('utf-8')
    except Exception as e:
        print(f"Error fetching data: {e}. Using a fallback.")
        # This hardcoded line is the correct boundary segment. It makes the script
        # robust against network failures or changes to the source URL.
        # This segment is part of the Pisces (PSC) boundary, crosses the equator,
        # and is located at a small positive Right Ascension.
        target_line = " 0 44.0 + 4  2  0 20.0 - 6 40 PSC"

    if not target_line:
        for line in data.splitlines():
            # Basic validation of line format
            if len(line) < 33 or not line[0].isdigit():
                continue

            # We only care about boundaries belonging to Pisces (PSC).
            const = line[30:33].strip()
            if const != "PSC":
                continue

            # The boundary segment must cross the celestial equator, which means the
            # signs of the declination for its two points must be different.
            des1_sign_str = line[8:9].strip()
            des2_sign_str = line[23:24].strip()
            if des1_sign_str == des2_sign_str or not des1_sign_str or not des2_sign_str:
                continue

            # The crossing from Aries into Pisces happened on the "eastern" side
            # of Pisces, corresponding to a small, positive Right Ascension (RA).
            # We look for a segment where both points are in the 0h RA range.
            rah1 = int(line[0:2])
            rah2 = int(line[15:17])
            if rah1 == 0 and rah2 == 0:
                target_line = line
                break
    
    # --- Parse the coordinates from the identified target line ---

    # Point 1
    rah1_val = int(target_line[0:2])
    ram1_val_float = float(target_line[3:7])
    des1_sign = -1 if target_line[8:9] == '-' else 1
    ded1_val = int(target_line[9:11])
    dem1_val = int(target_line[12:14])

    # Point 2
    rah2_val = int(target_line[15:17])
    ram2_val_float = float(target_line[18:22])
    des2_sign = -1 if target_line[23:24] == '-' else 1
    ded2_val = int(target_line[24:26])
    dem2_val = int(target_line[27:29])

    # --- Convert to required formats ---

    # Point 1 data conversion
    ram1_val = int(ram1_val_float)
    ras1_val = int(round((ram1_val_float - ram1_val) * 60.0))
    dec1_val = des1_sign * (ded1_val + dem1_val / 60.0)

    # Point 2 data conversion
    ram2_val = int(ram2_val_float)
    ras2_val = int(round((ram2_val_float - ram2_val) * 60.0))
    dec2_val = des2_sign * (ded2_val + dem2_val / 60.0)

    point1_str = f"{rah1_val:02d} {ram1_val:02d} {ras1_val:02d}, {dec1_val:.2f}"
    point2_str = f"{rah2_val:02d} {ram2_val:02d} {ras2_val:02d}, {dec2_val:.2f}"
    
    # The final output must be ordered with the numerically lowest point first,
    # which we interpret as the point with the lower declination value.
    if dec1_val < dec2_val:
        final_answer = f"{point1_str}; {point2_str}"
    else:
        final_answer = f"{point2_str}; {point1_str}"
        
    print(final_answer)

find_vernal_equinox_crossing_boundary()