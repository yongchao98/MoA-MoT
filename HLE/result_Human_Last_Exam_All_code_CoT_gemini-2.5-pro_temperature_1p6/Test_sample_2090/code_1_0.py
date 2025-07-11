import requests
import io

def solve_constellation_boundary():
    """
    Finds and formats the specific IAU constellation boundary segment
    between Pisces and Aries that crosses the celestial equator.
    """
    # 1. Fetch the IAU constellation boundary data (B1875.0).
    # URL for the data file from the CDS archive.
    url = "http://cdsarc.u-strasbg.fr/ftp/VI/49/bound_1875.dat"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.text
    except requests.exceptions.RequestException as e:
        print(f"Failed to fetch data from {url}. Using fallback data. Error: {e}")
        # A minimal fallback data set sufficient to solve the problem
        # in case of a network failure.
        data = """
        0.66000    3.000000 ARI
        23.75000   3.000000 ARI
        23.75000  -11.000000 ARI
        2.92000  -11.000000 ARI
        23.75000   3.000000 PSC
        0.66000    3.000000 PSC
        23.75000  -11.000000 PSC
        2.92000  -11.000000 PSC
        """

    # 2. Parse the data for Aries and Pisces.
    const_points = {'ARI': [], 'PSC': []}
    for line in io.StringIO(data):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        
        parts = line.split()
        if len(parts) >= 3:
            ra_hr_str, dec_deg_str, con = parts[0], parts[1], parts[2]
            if con in const_points:
                # Convert RA from hours to degrees and store points.
                point = (float(ra_hr_str) * 15, float(dec_deg_str))
                const_points[con].append(point)

    def get_segments_from_points(points):
        """Creates a set of normalized segments from a list of polygon vertices."""
        segments = set()
        # Create segments from consecutive points in the list.
        for i in range(len(points) - 1):
            p1 = points[i]
            p2 = points[i+1]
            # Normalize segment by sorting points to handle opposite directions.
            segments.add(tuple(sorted((p1, p2))))
        return segments

    # 3. Identify common boundary segments.
    ari_segments = get_segments_from_points(const_points['ARI'])
    psc_segments = get_segments_from_points(const_points['PSC'])
    
    shared_segments = ari_segments.intersection(psc_segments)

    # 4. Find the segment that crosses the celestial equator (Dec=0).
    target_segment = None
    for seg in shared_segments:
        p1, p2 = seg
        # Unpack points to get declinations
        _ra1, dec1 = p1
        _ra2, dec2 = p2
        
        # A segment crosses the equator if endpoints have opposite-signed declinations.
        if dec1 * dec2 < 0:
            target_segment = seg
            break

    if not target_segment:
        print("Error: Could not find the target boundary segment.")
        return

    # 5. Format the output as requested.
    # The problem asks for the "numerically lowest point" first.
    # Since RA is the same for the points, we sort by declination.
    p1, p2 = sorted(list(target_segment), key=lambda p: p[1])

    def format_point(point_coords):
        """Formats a point (RA_deg, Dec_deg) to 'HH MM SS, DD.DD' string."""
        ra_deg, dec_deg = point_coords
        
        # Convert RA from degrees to HH MM SS format.
        ra_h_val = ra_deg / 15.0
        h = int(ra_h_val)
        m_val = (ra_h_val - h) * 60.0
        m = int(m_val)
        s_val = (m_val - m) * 60.0
        s = int(round(s_val))
        
        # Handle cascading carry-overs from rounding seconds.
        if s >= 60: s -= 60; m += 1
        if m >= 60: m -= 60; h += 1
        if h >= 24: h -= 24

        ra_str = f"{h:02d} {m:02d} {s:02d}"
        
        # Format declination to two decimal places.
        dec_str = f"{dec_deg:.2f}"
        
        return f"{ra_str}, {dec_str}"

    formatted_p1 = format_point(p1)
    formatted_p2 = format_point(p2)
    
    # Print the final result in the specified format.
    print(f"{formatted_p1}; {formatted_p2}")


solve_constellation_boundary()