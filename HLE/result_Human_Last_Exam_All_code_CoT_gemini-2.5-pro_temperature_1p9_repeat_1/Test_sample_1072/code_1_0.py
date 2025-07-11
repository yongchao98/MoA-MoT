# Constants
surface_T = 25  # Surface temperature in °C
geo_gradient = 25  # Geothermal gradient in °C/km
zircon_Tc = 190  # Closure temperature for Zircon He in °C
apatite_Tc = 75   # Closure temperature for Apatite He in °C

# --- Sample 1: Zircon in Pluton ---
print("--- Sample 1 (Zircon) Age Calculation ---")
s1_start_depth = 15  # km
s1_start_time = 100  # Ma
s1_exhumation_rate = s1_start_depth / s1_start_time # km/Ma
print(f"Exhumation Rate: {s1_start_depth} km / {s1_start_time} Ma = {s1_exhumation_rate:.2f} km/Ma")

# Cooling rate in °C/Ma is Exhumation Rate * Geothermal Gradient
s1_cooling_rate_per_Ma = s1_exhumation_rate * geo_gradient
print(f"Cooling Rate: {s1_exhumation_rate:.2f} km/Ma * {geo_gradient} °C/km = {s1_cooling_rate_per_Ma:.2f} °C/Ma")

# Age is the time it takes to cool from Tc to surface T. Or more simply,
# find the time 't' when the sample was at depth D, where T(D) = Tc.
# T(t) = surface_T + (exhumation_rate * t) * geo_gradient
# Tc = surface_T + (rate * t * gradient)
# t = (Tc - surface_T) / (rate * gradient)
age_s1 = (zircon_Tc - surface_T) / (s1_exhumation_rate * geo_gradient)
print(f"Equation: Age = (Tc - Ts) / (Rate * Gradient)")
print(f"Calculation: Age = ({zircon_Tc} - {surface_T}) / ({s1_exhumation_rate:.2f} * {geo_gradient}) = {age_s1:.1f} Ma")
print(f"Sample 1 Age is approximately {age_s1:.1f} Ma.\n")


# --- Sample 2: Apatite in Sedimentary Rock ---
print("--- Sample 2 (Apatite) Age Calculation ---")
s2_start_T = 250 # °C
s2_start_time = 100 # Ma

# First, find the initial depth at 100 Ma
s2_start_depth = (s2_start_T - surface_T) / geo_gradient
print(f"Initial Depth at 100 Ma: ({s2_start_T} °C - {surface_T} °C) / {geo_gradient} °C/km = {s2_start_depth:.1f} km")
s2_exhumation_rate = s2_start_depth / s2_start_time
print(f"Exhumation Rate: {s2_start_depth:.1f} km / {s2_start_time} Ma = {s2_exhumation_rate:.2f} km/Ma")

# Now calculate the age using the same method as for sample 1
age_s2 = (apatite_Tc - surface_T) / (s2_exhumation_rate * geo_gradient)
print(f"Equation: Age = (Tc - Ts) / (Rate * Gradient)")
print(f"Calculation: Age = ({apatite_Tc} - {surface_T}) / ({s2_exhumation_rate:.2f} * {geo_gradient}) = {age_s2:.1f} Ma")
print(f"Sample 2 Age is approximately {age_s2:.1f} Ma.\n")


# --- Sample 3: Apatite in Rhyolite ---
print("--- Sample 3 (Apatite) Age ---")
age_s3 = 90 # Ma
print(f"Sample 3 erupted and cooled rapidly at 90 Ma.")
print(f"Sample 3 Age is approximately {age_s3} Ma.\n")

# --- Age Comparison ---
print("--- Age Comparison ---")
print(f"Sample 1: ~{age_s1:.0f} Ma")
print(f"Sample 2: ~{age_s2:.0f} Ma")
print(f"Sample 3: ~{age_s3:.0f} Ma")
print("\nConclusion: Sample 3 is the oldest, and Sample 2 is the youngest.")
print("This confirms that statement [H] is TRUE.")
