import math

# --- Given Parameters ---
num_levels = 6
largest_level_size_gb = 1.0  # GB
mem_buffer_size_kb = 1.0      # KB
insert_rate_bps = 16000.0     # bytes/s
page_size_bytes = 2500.0      # bytes

# --- Unit Conversions ---
# Convert all sizes to bytes for consistent calculations.
mem_buffer_size_bytes = mem_buffer_size_kb * 1024
largest_level_size_bytes = largest_level_size_gb * 1024**3

# --- Step 1: Calculate the Size Ratio (T) ---
# The number of levels on disk is num_levels - 1. The exponent is this value.
# Size(L_largest) = Size(L_memory) * T^(num_levels - 1)
num_disk_levels = num_levels - 1
size_ratio_T = (largest_level_size_bytes / mem_buffer_size_bytes)**(1.0 / num_disk_levels)

# --- Step 2: Calculate Total I/O Rate in Bytes/Second ---
# I/O Rate = (L0->L1 flush rate) + (Number of disk merges) * (Merge I/O rate)
# Flush Rate = R
# Merge Rate = R (read L_i) + R*T (read L_i+1) + R*T (write L_i+1) = R * (1 + 2T)
# Number of disk merges = num_disk_levels - 1 = 4
# Total Rate = R + 4 * R * (1 + 2T) = R * (1 + 4*(1+2T)) = R * (5 + 8T)
total_io_rate_bps = insert_rate_bps * (5 + 8 * size_ratio_T)

# --- Step 3: Convert to Page I/O Rate ---
page_io_rate = total_io_rate_bps / page_size_bytes

# --- Step 4: Output the results and the final equation ---
print("--- Calculation Steps ---")
print(f"1. Calculated Size Ratio (T): {size_ratio_T:.1f}")

# The question asks to output the numbers in the final equation.
# Final Equation: (Insert_Rate * (5 + 8 * T)) / Page_Size
print("\n2. Final Equation for Page I/O Rate:")
print(f"({insert_rate_bps:.1f} bytes/s * (5 + 8 * {size_ratio_T:.1f})) / {page_size_bytes:.1f} bytes/page")

# Calculate and print the final result
print(f"\n--- Result ---")
print(f"The minimum total page I/O rate is: {page_io_rate:.1f} pages/s")