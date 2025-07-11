import math

def titan_solver():
    """
    A program to simulate the Titan 6-bit computer's attempt to solve the physics problem.
    It demonstrates that the calculation is not feasible due to integer overflow.
    """

    # The maximum value for a 6-bit unsigned integer
    MAX_VAL = 63

    class TitanFrac:
        """
        Represents a number in Titan's fractional system: num/den * 10^exp.
        It checks for overflows during operations as per the architecture rules.
        """
        def __init__(self, num, den=1, exp=0):
            if den == 0:
                raise ValueError("Denominator cannot be zero.")
            # Simplification is part of the architecture to keep numbers small.
            common = math.gcd(int(num), int(den))
            self.num = int(num // common)
            self.den = int(den // common)
            self.exp = int(exp)

        def __mul__(self, other):
            """Simulates the MUL instruction, checking for overflow immediately."""
            print(f"Executing: ({self}) * ({other})")
            
            new_num = self.num * other.num
            new_den = self.den * other.den
            new_exp = self.exp + other.exp

            print(f"  - Intermediate product: {new_num}/{new_den} * 10^{new_exp}")

            if abs(new_num) > MAX_VAL or abs(new_den) > MAX_VAL:
                print(f"  - FAILURE: Overflow detected. The value '{new_num}' exceeds the 6-bit limit of {MAX_VAL}.")
                return None
            
            result = TitanFrac(new_num, new_den, new_exp)
            print(f"  - SUCCESS: Result is {result}")
            return result

        def __str__(self):
            return f"{self.num}/{self.den}e{self.exp}"

    print("--- Titan 6-bit Computer Feasibility Test ---")
    print("Goal: Calculate the mass of exoplanet Pandora.")
    print("Formula: Mass = Density * Volume = ρ * (4/3) * π * R³\n")

    # 1. Define constants in Titan's format
    print("Step 1: Defining constants.")
    try:
        # Density (ρ): 1.2 metric tons/m³ = 1200 kg/m³ -> 12/1 * 10^2
        rho = TitanFrac(12, 1, 2)
        print(f"  - Density (ρ): {rho}")

        # Radius (R): 2000 km = 2,000,000 m -> 2/1 * 10^6
        R_planet = TitanFrac(2, 1, 6)
        print(f"  - Radius (R): {R_planet}")
        
        # Constant 4/3
        four_thirds = TitanFrac(4, 3)
        print(f"  - Constant: {four_thirds}")
        
        # Using a simple approximation for Pi (π ≈ 3)
        pi_approx = TitanFrac(3, 1)
        print(f"  - Pi approximation (π): {pi_approx}\n")

    except ValueError as e:
        print(f"Initialization Error: {e}")
        print("<<<N0>>>")
        return

    # 2. Calculate Volume components
    print("Step 2: Calculating Volume components (V = 4/3 * π * R³).")
    
    # R³ = R * R * R
    R_squared = R_planet * R_planet
    R_cubed = R_squared * R_planet
    if R_cubed is None: # This case won't be hit, but good practice
        print("\nFATAL ERROR: Overflow during radius calculation.")
        print("<<<N0>>>")
        return
    print(f"-> Calculated R³ = {R_cubed}\n")

    # 3. Calculate Mass step-by-step
    print("Step 3: Calculating Mass = (ρ * 4/3 * π) * R³.")
    
    # Mass = ρ * (4/3)
    mass_term = rho * four_thirds
    if mass_term is None: # This case won't be hit, but good practice
        print("\nFATAL ERROR: Overflow during mass calculation.")
        print("<<<N0>>>")
        return

    # Mass = (ρ * 4/3) * π
    mass_term = mass_term * pi_approx
    if mass_term is None: # This case won't be hit, but good practice
        print("\nFATAL ERROR: Overflow during mass calculation.")
        print("<<<N0>>>")
        return
        
    # Mass = (ρ * 4/3 * π) * R³
    # This is the step that will fail.
    # Current mass_term is 48/1e2. R_cubed is 8/1e18.
    # The multiplication will be 48 * 8 = 384, which is > 63.
    final_mass = mass_term * R_cubed

    if final_mass is None:
        print("\n--- CONCLUSION ---")
        print("The calculation is NOT POSSIBLE on the Titan architecture.")
        print("The multiplication of core physical constants (like density, radius, pi) results in numerators that exceed the 6-bit limit.")
        print("Since the planet's mass cannot be computed, subsequent calculations for the Schwarzschild radius and gravitational force are also impossible.")
        print("\nFinal Answer:")
        print("<<<N0>>>")
    else:
        # This code path is unreachable
        print(f"Final Mass: {final_mass}")
        print("<<<Y?>>>")

titan_solver()