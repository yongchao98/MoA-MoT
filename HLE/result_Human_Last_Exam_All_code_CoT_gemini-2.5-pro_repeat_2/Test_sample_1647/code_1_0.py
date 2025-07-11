import math

# Step 1: Define constants and problem values.
# The distance 'd' from the event horizon is given.
d_km = 13

# The Schwarzschild Radius 'Rs' must be calculated from Pandora's orbital data
# and then approximated for the Bagua 'frac' type.
# Rs = (8 * pi^2 * a^3) / (T^2 * c^2)
# a = 1e11 m, T = 800 days = 69,120,000 s, c = 299792458 m/s
# Rs ≈ 184 meters = 0.184 km.
# To fit this into a Bagua 'frac' (n/d*10^e) where 'n' is a signed char (-32 to 31),
# we can use the literal 18/10e-1, which evaluates to 0.18.
Rs_km_approx = 0.18

# Step 2: Calculate the time dilation factor 'f' and round it.
f_value = math.sqrt(d_km / (Rs_km_approx + d_km))
f_rounded = round(f_value, 3)

# Step 3: Calculate the memory usage 'z' for the most memory-efficient Bagua C program.
# The program needs variables for the input 'd' and for storing results.
# - 'd_km': Can be stored in a 'signed char' (13 is in the range -32 to 31). Size = 2 trits.
# - 'f_squared': Stores the intermediate fraction d/(Rs+d). Type 'frac'. Size = 8 trits.
# - 'f': Stores the final result after the sqrt approximation. Type 'frac'. Size = 8 trits.
# The Rs value is a constant literal in the code, not a variable.
mem_d_km = 2
mem_f_squared = 8
mem_f = 8
z_trits = mem_d_km + mem_f_squared + mem_f

# Step 4: Print the calculation steps and the final answer.
print("Time Dilation Factor (f) Calculation:")
print(f"f = sqrt(d / (Rs + d))")
print(f"f = sqrt({d_km} / ({Rs_km_approx} + {d_km}))")
print(f"f ≈ {f_rounded}\n")

print("Memory Usage (z) Calculation:")
print(f"z = size_of(d_km) + size_of(f_squared) + size_of(f)")
print(f"z = {mem_d_km} trits (char) + {mem_f_squared} trits (frac) + {mem_f} trits (frac)")
print(f"z = {z_trits} trits\n")

# Step 5: Output the final answer in the required format f:z.
final_answer = f"{f_rounded}:{z_trits}"
print(f"<<<{final_answer}>>>")