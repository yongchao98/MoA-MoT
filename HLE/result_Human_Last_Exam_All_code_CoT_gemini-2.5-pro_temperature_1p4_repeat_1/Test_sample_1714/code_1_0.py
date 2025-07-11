import math

# This Python script simulates the calculation on the Wuxing architecture
# to find the dark matter percentage of the Pandora galaxy.

def calculate_dark_matter_percentage():
    """
    Performs the dark matter calculation based on the problem's parameters
    and the constraints of the Wuxing computing environment.
    """
    
    # --- Part 1: Define Constants in Wuxing frac Representation ---
    # Using astronomical units (M_sun, kpc, km/s) to keep numbers manageable.
    # frac = (numerator, denominator, exponent)
    
    # v = 200 km/s -> (2/1) * 10^2
    v = (2, 1, 2)
    # r = 10 kpc -> (1/1) * 10^1
    r = (1, 1, 1)
    # M_lum = 6e9 M_sun -> (6/1) * 10^9
    M_lum = (6, 1, 9)
    # G_astro ≈ 4.3e-6 -> (43/10) * 10^-6
    G = (43, 10, -6)
    
    # --- Part 2: Simulate the Wuxing `frac` Calculation ---
    
    # M_total = (v^2 * r) / G
    # Step-by-step to show the logic:
    # v^2 = {4, 1, 4}
    # v^2 * r = {4, 1, 5}
    # M_total = {4, 1, 5} / {43, 10, -6} => {4*10, 1*43, 5-(-6)} = {40, 43, 11}
    M_total = (40, 43, 11)

    # ratio = M_luminous / M_total
    # Raw calculation: {6,1,9} / {40,43,11} => {6*43, 1*40, 9-11} = {258, 40, -2}
    # Simplified: {129, 20, -2}. Numerator 129 overflows signed char.
    # Here we simulate the "overflow protection" by approximating the value.
    # Value is (129/20) * 10^-2 = 0.0645.
    # We approximate to 0.064, which is 64/1000 = 8/125.
    # This becomes the storable frac {8, 125, 0}.
    ratio = (8, 125, 0)
    
    # percentage = (1 - ratio) * 100
    # 1 - {8, 125, 0} = (125-8)/125 = {117, 125, 0}
    one_minus_ratio = (117, 125, 0)
    # Multiply by 100 ({1, 1, 2}): {117, 125, 0} -> {117, 125, 2}
    percentage_frac = (117, 125, 2)

    # Convert final frac to a floating point number for the answer
    p = (percentage_frac[0] / percentage_frac[1]) * (10**percentage_frac[2])

    # --- Part 3: Calculate Memory Usage (z) ---
    # The hypothetical C program would need:
    # 7 frac variables (v, r, G, M_lum, M_total, ratio, percentage)
    # 3 int variables for printing
    # Memory: frac = 6D, int = 5D
    z = (7 * 6) + (3 * 5)
    
    # --- Part 4: Output the Final Answer ---
    
    print("Final Equation:")
    print(f"(1 - (M_luminous / M_total)) * 100")
    print(f"= (1 - ({M_lum[0]}e{M_lum[2]} / ({M_total[0]}/{M_total[1]})e{M_total[2]})) * 100")
    print(f"= (1 - {ratio[0]}/{ratio[1]}) * 100")
    print(f"= {p:.1f}%\n")

    print("Final Answer (p:z):")
    print(f"{p:.1f}:{z}")

calculate_dark_matter_percentage()