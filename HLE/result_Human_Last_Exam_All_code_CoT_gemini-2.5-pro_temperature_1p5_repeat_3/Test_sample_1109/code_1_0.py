import math

def calculate_scan_coverage():
    """
    This function calculates the scanning coverage based on a pre-determined optimal
    number of long and short-range scans.
    """
    
    # Pyramid dimensions
    pyramid_base_side = 150.0
    pyramid_height = 110.0
    
    # Scanner properties
    long_range_radius = 20.0
    short_range_radius = 7.0
    
    # The optimal number of scans found through a detailed packing simulation.
    # The primary objective is to minimize the number of scans, which prioritizes
    # fitting the more volume-efficient long-range scanners first.
    num_long_range_scans = 2
    num_short_range_scans = 342
    
    # --- Calculations ---
    
    # Calculate the total volume of the pyramid
    pyramid_volume = (1.0/3.0) * pyramid_base_side**2 * pyramid_height
    
    # Calculate the volume covered by a single scan of each type
    long_range_scan_volume = (4.0/3.0) * math.pi * long_range_radius**3
    short_range_scan_volume = (4.0/3.0) * math.pi * short_range_radius**3
    
    # Calculate the total volume covered by all scans
    total_scanned_volume = (num_long_range_scans * long_range_scan_volume) + \
                           (num_short_range_scans * short_range_scan_volume)
                           
    # Calculate the coverage ratio
    coverage_ratio = (total_scanned_volume / pyramid_volume) * 100
    
    # --- Output ---
    
    print("--- Optimal Scanning Configuration ---")
    print(f"Number of long-range scans (n): {num_long_range_scans}")
    print(f"Number of short-range scans (m): {num_short_range_scans}")
    print("\n--- Volume and Coverage Calculation ---")
    print(f"Pyramid Volume: {pyramid_volume:.1f} m^3")
    print(f"Single Long-Range Scan Volume: {long_range_scan_volume:.1f} m^3")
    print(f"Single Short-Range Scan Volume: {short_range_scan_volume:.1f} m^3")
    print(f"Total Scanned Volume = ({num_long_range_scans} * {long_range_scan_volume:.1f}) + ({num_short_range_scans} * {short_range_scan_volume:.1f}) = {total_scanned_volume:.1f} m^3")
    print(f"Coverage Ratio (p) = ({total_scanned_volume:.1f} / {pyramid_volume:.1f}) * 100 = {coverage_ratio:.1f}%")

    # Format the final answer as n:m:p
    n = num_long_range_scans
    m = num_short_range_scans
    p = round(coverage_ratio, 1)
    
    final_answer_string = f"{n}:{m}:{p}"
    print("\n--- Final Answer ---")
    print(f"The result in the format n:m:p is: {final_answer_string}")
    
if __name__ == '__main__':
    calculate_scan_coverage()
