import math

# Step 1: Define constants and parameters
L = 6  # Number of levels
S_L_gb = 1  # Largest level size in GB
S_mem_kb = 1  # Memory buffer size in KB
R_insert = 16000  # Insert rate in bytes/s
P = 2500  # Page size in bytes

# Convert sizes to bytes for consistent units
S_L = S_L_gb * (1024**3)  # 1 GB in bytes
S_mem = S_mem_kb * 1024    # 1 KB in bytes

print("--- Input Parameters ---")
print(f"Number of levels (L): {L}")
print(f"Largest level size (S_L): {S_L_gb} GB ({S_L} bytes)")
print(f"Memory buffer size (S_mem): {S_mem_kb} KB ({S_mem} bytes)")
print(f"Insert rate (R_insert): {R_insert} bytes/s")
print(f"Page size (P): {P} bytes")
print("-" * 26)

# Step 2: Calculate the size ratio (T)
# S_L = S_mem * T^L  =>  T = (S_L / S_mem)^(1/L)
T = (S_L / S_mem)**(1 / L)

print("\n--- Calculation Steps ---")
print(f"1. Calculate Size Ratio (T):")
print(f"   T = (S_L / S_mem)^(1/L)")
print(f"   T = ({S_L} / {S_mem})^(1/{L})")
print(f"   T = {T:.4f}")

# Step 3: Calculate the total I/O rate in bytes/s
# Total Write Rate = L * R_insert
# Total Read Rate = (L-1) * (1 + T) * R_insert
# Total I/O Rate = Write Rate + Read Rate = R_insert * (L + (L-1)*(1+T))
total_io_bytes_per_sec = R_insert * (L + (L - 1) * (1 + T))

print(f"\n2. Calculate Total I/O Rate in bytes/s:")
print(f"   Formula: R_insert * (L + (L-1) * (1 + T))")
print(f"   Calculation: {R_insert} * ({L} + ({L-1}) * (1 + {T:.4f}))")
print(f"   Total I/O Rate = {total_io_bytes_per_sec:.4f} bytes/s")

# Step 4: Convert to page I/O rate
# Page I/O Rate = Total I/O Rate (bytes/s) / Page Size (bytes)
page_io_rate = total_io_bytes_per_sec / P

print(f"\n3. Calculate Total Page I/O Rate:")
print(f"   Formula: Total I/O Rate (bytes/s) / Page Size")
print(f"   Calculation: {total_io_bytes_per_sec:.4f} / {P}")
print("-" * 26)
print(f"\nFinal Answer: The minimum total page I/O rate is {page_io_rate:.4f} pages/s.")
