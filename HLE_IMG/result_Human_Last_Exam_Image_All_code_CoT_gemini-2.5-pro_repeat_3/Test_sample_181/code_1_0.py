import math

# Step 1: Define physical constants
# Average distance from Earth to Sun in meters
R_es = 1.496e11
# Earth's orbital period in seconds (365.25 days)
T_e = 365.25 * 24 * 60 * 60

# Average distance from Moon to Earth in meters
R_me = 3.844e8
# Moon's sidereal orbital period in seconds (27.3 days)
T_m = 27.3 * 24 * 60 * 60

# Step 2: Calculate Earth's orbital speed around the Sun
v_e = (2 * math.pi * R_es) / T_e

# Step 3: Calculate Moon's orbital speed around the Earth
v_m = (2 * math.pi * R_me) / T_m

# Step 4 & 5: Compare the speeds and print the results
print("Analyzing the Moon's orbit around the Sun:")
print("-" * 40)
print(f"Earth's orbital speed around the Sun (v_e):")
print(f"  v_e = (2 * pi * R_es) / T_e")
print(f"  v_e = (2 * {math.pi:.4f} * {R_es:.3e} m) / ({T_e:.3e} s)")
print(f"  v_e = {v_e:,.2f} m/s (or {v_e/1000:,.2f} km/s)")
print("-" * 40)
print(f"Moon's orbital speed around the Earth (v_m):")
print(f"  v_m = (2 * pi * R_me) / T_m")
print(f"  v_m = (2 * {math.pi:.4f} * {R_me:.3e} m) / ({T_m:.3e} s)")
print(f"  v_m = {v_m:,.2f} m/s (or {v_m/1000:,.2f} km/s)")
print("-" * 40)
print("Conclusion:")
if v_e > v_m:
    print(f"Since the Earth's speed ({v_e/1000:,.2f} km/s) is much greater than the Moon's speed ({v_m/1000:,.2f} km/s),")
    print("the Moon is always moving forward in its orbit around the Sun.")
    print("Its path is a wave, but it never loops backward.")
    print("Therefore, the orbit looks most like option C.")
else:
    print(f"Since the Moon's speed ({v_m/1000:,.2f} km/s) is greater than the Earth's speed ({v_e/1000:,.2f} km/s),")
    print("the Moon's path would have loops or cusps.")
    print("This would look like option D.")
