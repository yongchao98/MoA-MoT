import math

# --- User-defined parameters ---

# Time instant up to which convergence is required
T = 20.0 

# Propagation speed of the wave
c = 50.0

# Right boundary of the left subdomain Omega_1
b = 12.0

# Left boundary of the right subdomain Omega_2
a = 10.0

# --- Calculation ---

# Calculate the size of the overlap
# M = b - a
M = b - a

# Check for valid overlap
if M <= 0:
    print("Error: The domains do not overlap (b must be greater than a).")
else:
    # Calculate the ratio Tc/M
    ratio = (T * c) / M
    
    # The number of updates required is the ceiling of this ratio
    num_updates = math.ceil(ratio)
    
    # The iteration counter starts at 0 and is incremented after each update.
    # So, after N updates, the counter shows N-1.
    iterations_needed = num_updates - 1

    # --- Output the result ---
    print(f"Given parameters:")
    print(f"Time Horizon (T) = {T}")
    print(f"Wave Speed (c) = {c}")
    print(f"Overlap Start (a) = {a}")
    print(f"Overlap End (b) = {b}")
    print(f"Overlap Size (M = b - a) = {M}\n")

    print(f"Calculation:")
    print(f"Number of iterations = ceil(T * c / M) - 1")
    print(f"                     = ceil({T} * {c} / {M}) - 1")
    print(f"                     = ceil({ratio}) - 1")
    print(f"                     = {num_updates} - 1")
    print(f"                     = {iterations_needed}\n")

    print(f"Final answer: The method needs to perform iterations 0, 1, ..., {iterations_needed}.")
    print(f"The final value on the iteration counter will be: {iterations_needed}")
