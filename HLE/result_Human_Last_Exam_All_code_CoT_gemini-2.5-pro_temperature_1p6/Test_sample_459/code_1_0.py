import math

# --- Given Parameters ---
num_levels = 6
size_lmax_gb = 1.0
size_mem_kb = 1.0
insert_rate_bps = 16000.0  # bytes per second
page_size_bytes = 2500.0

# --- Step 1: Convert all sizes to bytes ---
size_lmax_bytes = size_lmax_gb * (1024**3)
size_mem_bytes = size_mem_kb * 1024

# --- Step 2: Calculate the Size Ratio (T) ---
# Size_Lmax = T^num_levels * Size_mem
# T = (Size_Lmax / Size_mem) ^ (1 / num_levels)
size_ratio_base = size_lmax_bytes / size_mem_bytes
T = size_ratio_base ** (1.0 / num_levels)

# --- Step 3: Calculate Total Byte I/O Rate ---
# This is for leveled compaction, which gives the minimum I/O.
# For each unit of data moved from L(i) to L(i+1), we have:
# - 1 read from L(i)
# - T reads from L(i+1)
# - T+1 writes to L(i+1)
# Total I/O multiplier per compaction step = 1 + T + (T+1) = 2T + 2
compaction_io_multiplier = 2 * T + 2
num_compactions = num_levels - 1

# Total I/O Rate = (Flush I/O) + (Compaction I/O)
# Total I/O Rate = insert_rate + num_compactions * insert_rate * (2T+2)
total_byte_io_rate = insert_rate_bps * (1 + num_compactions * compaction_io_multiplier)

# --- Step 4: Calculate Total Page I/O Rate ---
total_page_io_rate = total_byte_io_rate / page_size_bytes

# --- Final Output ---
print("--- Calculation Steps ---")
print(f"1. Number of Levels (L): {num_levels}")
print(f"2. Size Ratio (T): {T:.4f}")
print(f"3. Insert Rate (bytes/s): {insert_rate_bps}")
print(f"4. Page Size (bytes): {page_size_bytes}")
print("\n--- Final Equation ---")
# Remember to output each number in the final equation!
equation_str = (
    f"Total Page I/O Rate = "
    f"({insert_rate_bps:.0f} bytes/s * (1 + ({num_levels} - 1) * (2 * {T:.4f} + 2))) / {page_size_bytes:.0f} bytes/page"
)
print(equation_str)
print("\n--- Result ---")
print(f"The minimum total page I/O rate is: {total_page_io_rate:.2f} pages/s")

# Final answer in the required format
#<<<715.48>>>