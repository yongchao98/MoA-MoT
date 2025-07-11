import math

def calculate_lithiation_time():
    """
    Calculates the lithiation time for a larger, quarter-sphere LFP particle
    based on the lithiation time of a smaller, spherical particle.
    """
    # --- Case 1: Initial Particle ---
    # A full sphere
    d1_nm = 50.0  # Diameter in nanometers
    t1_s = 20.0   # Lithiation time in seconds

    # --- Case 2: New Particle ---
    # A quarter of a sphere
    d2_nm = 250.0 # Diameter in nanometers

    # Step 1: Calculate the radius for both particles.
    # The characteristic diffusion length for both a sphere and a quarter-sphere is its radius.
    r1_nm = d1_nm / 2
    r2_nm = d2_nm / 2

    # Step 2: Apply the diffusion time relationship.
    # The time for a diffusion-limited process (t) is proportional to the square of the
    # characteristic length (r), i.e., t ∝ r².
    # Therefore, we can set up the ratio: t2 / t1 = (r2 / r1)²
    # We solve for t2: t2 = t1 * (r2 / r1)²
    t2_s = t1_s * (r2_nm / r1_nm)**2

    # Step 3: Print the results and the equation with the final numbers.
    print("--- Initial Particle Parameters ---")
    print(f"Shape: Sphere")
    print(f"Diameter (d₁): {d1_nm} nm")
    print(f"Radius (r₁): {r1_nm} nm")
    print(f"Lithiation Time (t₁): {t1_s} s\n")

    print("--- New Particle Parameters ---")
    print(f"Shape: Quarter-Sphere")
    print(f"Diameter (d₂): {d2_nm} nm")
    print(f"Radius (r₂): {r2_nm} nm\n")

    print("--- Calculation ---")
    print("The lithiation time (t) is proportional to the square of the radius (r).")
    print("The formula is: t₂ = t₁ * (r₂ / r₁)²\n")
    print("Plugging in the numbers:")
    # The final equation with numbers as requested
    print(f"t₂ = {t1_s} s * ({r2_nm} nm / {r1_nm} nm)²")
    ratio = r2_nm / r1_nm
    print(f"t₂ = {t1_s} s * ({ratio})²")
    print(f"t₂ = {t1_s} s * {ratio**2}")
    print(f"t₂ = {t2_s} s\n")
    print(f"The calculated lithiation time for the quarter-sphere particle is {t2_s} seconds.")
    return t2_s

# Execute the function and get the final answer.
final_time = calculate_lithiation_time()
# The prompt requires the answer in a specific format at the end.
# print(f"\n<<< {final_time} >>>")
# But let's just output the final value as a number.
final_answer = int(final_time)
# The final result is in the last line of output.
print(f"<<<{final_answer}>>>")
