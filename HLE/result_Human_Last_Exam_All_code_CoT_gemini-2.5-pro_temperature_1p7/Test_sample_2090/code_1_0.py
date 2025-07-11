import math

def format_boundary_points():
    """
    This function identifies the IAU constellation boundary between Pisces and Aries
    that crosses the celestial equator, formats the coordinates of its endpoints,
    and prints them in the specified format.
    """

    # The B1875.0 coordinates for the boundary segment vertices
    # Point 1 (lower declination)
    p1_ra_h_raw = 23
    p1_ra_m_raw = 59.2
    p1_dec_d_raw = -6
    p1_dec_m_raw = 11

    # Point 2 (higher declination)
    p2_ra_h_raw = 23
    p2_ra_m_raw = 59.2
    p2_dec_d_raw = 5
    p2_dec_m_raw = 35

    # --- Convert Point 1 ---
    # Convert RA from hours and decimal minutes to hours, minutes, seconds
    p1_ra_m = int(p1_ra_m_raw)
    p1_ra_s = round((p1_ra_m_raw - p1_ra_m) * 60)

    # Convert Dec from degrees and arcminutes to decimal degrees
    # The sign is determined by the degrees component
    p1_dec_decimal = p1_dec_d_raw - (p1_dec_m_raw / 60.0)

    # --- Convert Point 2 ---
    # RA is the same for both points
    p2_ra_m = int(p2_ra_m_raw)
    p2_ra_s = round((p2_ra_m_raw - p2_ra_m) * 60)

    # Convert Dec from degrees and arcminutes to decimal degrees
    p2_dec_decimal = p2_dec_d_raw + (p2_dec_m_raw / 60.0)
    
    # Per the instructions, outputting each number in the final line.
    # The components for the first point (lowest declination):
    p1_ra_h_final = p1_ra_h_raw
    p1_ra_m_final = p1_ra_m
    p1_ra_s_final = p1_ra_s
    p1_dec_final = p1_dec_decimal

    # The components for the second point:
    p2_ra_h_final = p2_ra_h_raw
    p2_ra_m_final = p2_ra_m
    p2_ra_s_final = p2_ra_s
    p2_dec_final = p2_dec_decimal

    # Print the components first, as requested by the prompt.
    print(f"Point 1 (Lowest Declination):")
    print(f"RA Hour: {p1_ra_h_final}")
    print(f"RA Minute: {p1_ra_m_final}")
    print(f"RA Second: {p1_ra_s_final}")
    print(f"Declination: {p1_dec_final:.2f}")

    print(f"\nPoint 2 (Highest Declination):")
    print(f"RA Hour: {p2_ra_h_final}")
    print(f"RA Minute: {p2_ra_m_final}")
    print(f"RA Second: {p2_ra_s_final}")
    print(f"Declination: {p2_dec_final:.2f}")

    # Format the final strings for each point
    p1_str = f"{p1_ra_h_final} {p1_ra_m_final:02d} {p1_ra_s_final:02d}, {p1_dec_final:+06.2f}"
    p2_str = f"{p2_ra_h_final} {p2_ra_m_final:02d} {p2_ra_s_final:02d}, {p2_dec_final:+06.2f}"
    
    # Combine them with a semicolon for the final answer
    final_answer = f"{p1_str}; {p2_str}"
    
    print("\nFinal Formatted Answer:")
    print(final_answer)

format_boundary_points()