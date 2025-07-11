import math

# Given parameters
num_disk_levels = 6
largest_level_size_gb = 1.0
memory_buffer_kb = 1.0
insert_rate_bps = 16000.0
page_size_bytes = 2500.0

# --- Calculations ---

# Convert all sizes to bytes for consistency
largest_level_size_bytes = largest_level_size_gb * (1024**3)
memory_buffer_bytes = memory_buffer_kb * 1024

# Step 1: Calculate the size ratio (T)
# With 1 memory level (L0) and 6 disk levels (L1-L6), there are 6 size steps.
# Size(L6) = T^6 * Size(L0)
num_size_steps = num_disk_levels
t_pow_n = largest_level_size_bytes / memory_buffer_bytes
size_ratio_T = t_pow_n**(1.0/num_size_steps)

# Step 2: Calculate the total I/O rate in bytes/second
# Model: Total I/O = Initial Write + Compaction I/O
# Initial Write Rate (L0->L1) = insert_rate_bps
# Compaction I/O: There are (num_disk_levels - 1) = 5 compactions.
# For each compaction, read amplification is T and write amplification is T.
# Total I/O Rate = Insert_Rate * (1 + (Num_levels - 1) * (T_read_amp + T_write_amp))
# Total I/O Rate = Insert_Rate * (1 + 5 * (T + T)) = Insert_Rate * (1 + 10 * T)
num_compactions = num_disk_levels - 1
total_io_rate_bps = insert_rate_bps * (1 + num_compactions * 2 * size_ratio_T)

# Step 3: Convert the I/O rate from bytes/s to pages/s
total_page_io_rate = total_io_rate_bps / page_size_bytes

# --- Output the explanation and final equation ---

print("### LSM Tree I/O Calculation ###")
print(f"Given Parameters:")
print(f"  - Disk Levels: {num_disk_levels}")
print(f"  - Memory Buffer (L0): {memory_buffer_kb} KB")
print(f"  - Largest Level (L6): {largest_level_size_gb} GB")
print(f"  - Insert Rate: {insert_rate_bps} bytes/s")
print(f"  - Page Size: {page_size_bytes} bytes")
print("-" * 30)

print("Step 1: Calculate the Size Ratio (T)")
print(f"The size ratio 'T' is calculated across {num_size_steps} levels (from L0 to L6).")
print(f"T = (Size(L6) / Size(L0)) ^ (1/{num_size_steps})")
print(f"T = ({largest_level_size_bytes:.0f} / {memory_buffer_bytes:.0f}) ^ (1/{num_size_steps}) = {t_pow_n:.0f} ^ {1.0/num_size_steps:.3f}")
print(f"T = {size_ratio_T:.4f}")
print("-" * 30)

print("Step 2: Calculate Total I/O Rate in Bytes/s")
print("Total I/O = (Initial Flush I/O) + (Compaction I/O)")
print(f"The model for minimum I/O (leveled compaction) is:")
print("  Total Rate = Insert Rate * (1 [for L0->L1 write] + Num_Compactions * (Read_Amp + Write_Amp))")
print(f"  Read/Write amplification per compaction level is T = {size_ratio_T:.4f}")
print(f"  Number of compactions between disk levels = {num_compactions}")
print(f"Total I/O Rate (bytes/s) = {insert_rate_bps:.0f} * (1 + {num_compactions} * ({size_ratio_T:.4f} + {size_ratio_T:.4f}))")
print(f"Total I/O Rate (bytes/s) = {insert_rate_bps:.0f} * (1 + {num_compactions*2} * {size_ratio_T:.4f})")
print(f"Total I/O Rate (bytes/s) = {total_io_rate_bps:.2f}")
print("-" * 30)

print("Step 3: Convert to Page I/O Rate")
print("Page I/O Rate = Total I/O Rate (bytes/s) / Page Size (bytes)")
print(f"Final Equation: {total_io_rate_bps:.2f} / {page_size_bytes:.0f}")
print(f"Minimum Total Page I/O Rate = {total_page_io_rate:.2f} pages/s")
