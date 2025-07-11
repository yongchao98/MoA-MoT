import urllib.request
import sys

def get_constellation_data():
    """
    Fetches and parses the B1875.0 IAU constellation boundary data.
    The data is parsed into a dictionary where keys are constellation
    abbreviations and values are lists of (RA_hours, Dec_degrees) tuples.
    """
    # URL for the official B1875.0 boundary data from CDS VizieR
    url = "http://cdsarc.u-strasbg.fr/ftp/VI/49/bounds.dat"
    try:
        response = urllib.request.urlopen(url)
        data = response.read().decode('utf-8')
    except Exception as e:
        print(f"Error fetching data: {e}", file=sys.stderr)
        return None

    lines = data.strip().split('\n')
    
    boundaries = {}
    current_con = ""
    for line in lines:
        line = line.strip()
        # Skip empty or header lines
        if not line or "---" in line or "RAh" in line:
            continue
        
        parts = line.split()
        
        # A single 3-letter part indicates a new constellation
        if len(parts) == 1 and len(parts[0]) == 3 and parts[0].isalpha():
            current_con = parts[0]
            if current_con not in boundaries:
                boundaries[current_con] = []
        # A line with coordinates
        elif len(parts) >= 3 and current_con:
            try:
                rah = int(parts[0])
                ram = int(parts[1])
                ded = int(parts[2])
                
                # Convert RA to decimal hours
                ra_h = rah + ram / 60.0
                # Declination is already in integer degrees
                dec_d = float(ded)
                
                boundaries[current_con].append((ra_h, dec_d))
            except (ValueError, IndexError):
                # Ignore malformed lines
                continue
                
    return boundaries

def find_boundary_marker():
    """
    Identifies the Aries/Pisces boundary segment that crosses the celestial equator.
    """
    all_boundaries = get_constellation_data()
    if not all_boundaries:
        print("Failed to retrieve boundary data.", file=sys.stderr)
        return

    ari_points = all_boundaries.get("ARI")
    psc_points = all_boundaries.get("PSC")

    if not ari_points or not psc_points:
        print("Could not find boundary data for Aries or Pisces.", file=sys.stderr)
        return

    # Create a set of directed segments for Aries for efficient lookup
    ari_segments = set()
    for i in range(len(ari_points)):
        p1 = ari_points[i]
        p2 = ari_points[(i + 1) % len(ari_points)]
        ari_segments.add((p1, p2))

    # Find the shared segment that crosses the celestial equator (Dec=0)
    target_segment = None
    for i in range(len(psc_points)):
        p1 = psc_points[i]
        p2 = psc_points[(i + 1) % len(psc_points)]
        
        # A segment is shared if the reverse segment exists in the other list
        if (p2, p1) in ari_segments:
            # Check if this segment crosses the equator (Dec=0)
            # This occurs if the declinations of the endpoints have opposite signs.
            if (p1[1] > 0 and p2[1] < 0) or (p1[1] < 0 and p2[1] > 0):
                target_segment = (p1, p2)
                break
    
    if not target_segment:
        print("Could not identify the specific boundary marker.", file=sys.stderr)
        return
        
    # Sort points by declination to ensure the lowest appears first
    p1, p2 = sorted(target_segment, key=lambda p: p[1])
    
    def format_point(point):
        """Formats a point into the specified string format."""
        ra_h, dec_d = point
        
        # Convert RA from decimal hours to hh mm ss
        h = int(ra_h)
        m_rem = (ra_h - h) * 60
        m = int(m_rem)
        s_rem = (m_rem - m) * 60
        s = int(round(s_rem))
        
        # Handle cases where rounding seconds increments minutes/hours
        if s == 60:
            s = 0
            m += 1
        if m == 60:
            m = 0
            h += 1

        ra_str = f"{h:02d} {m:02d} {s:02d}"
        
        # Format declination to two decimal places
        dec_str = f"{dec_d:.2f}"
        
        return f"{ra_str}, {dec_str}"

    point1_str = format_point(p1)
    point2_str = format_point(p2)
    
    print(f"{point1_str}; {point2_str}")


if __name__ == '__main__':
    find_boundary_marker()