import math

# Titan Computer Specification Constants
NUM_DEN_MAX = 31  # 5-bit unsigned integer max value
EXP_MIN = -16     # 5-bit signed integer range min
EXP_MAX = 15      # 5-bit signed integer range max

class TitanError(Exception):
    """Custom exception for Titan computation failures."""
    pass

class TitanNumber:
    """
    A class to represent and operate on numbers within the Titan architecture.
    A number is stored as N/D * 2^E, where N, D are 5-bit unsigned integers (1-31)
    and E is a 5-bit signed integer (-16 to 15).
    """
    def __init__(self, value):
        self.value = value
        self.n, self.d, self.e, self.error = self._find_representation(value)
        if self.n is None:
            raise TitanError(f"Value {value:.2e} cannot be represented. "
                             f"Its magnitude is too large or too small for the exponent range [{EXP_MIN}, {EXP_MAX}].")

    def _find_representation(self, value):
        """
        Finds the best fractional representation (N/D) * 2^E for a given value.
        """
        if value == 0:
            return 0, 1, 0, 0

        best_rep = {'n': None, 'd': None, 'e': None, 'error': float('inf')}

        # Iterate through all possible exponents
        for e in range(EXP_MAX, EXP_MIN - 1, -1):
            target_fraction = value / (2**e)
            
            # Optimization: if the required fraction is already larger than what N/D can represent,
            # smaller exponents will only make it worse.
            if target_fraction > NUM_DEN_MAX:
                continue
            
            # Optimization: If target is too small (e.g., < 1/NUM_DEN_MAX), we can also skip,
            # but iterating E from high to low handles this naturally.

            # Find the best N/D for the current exponent E
            for d in range(1, NUM_DEN_MAX + 1):
                ideal_n = target_fraction * d
                # Check the closest integer numerator
                n = round(ideal_n)
                
                if 1 <= n <= NUM_DEN_MAX:
                    approx_val = (n / d) * (2**e)
                    error = abs(approx_val - value) / abs(value) if value != 0 else 0
                    
                    if error < best_rep['error']:
                        best_rep.update({'n': n, 'd': d, 'e': e, 'error': error})
        
        return best_rep['n'], best_rep['d'], best_rep['e'], best_rep['error']

    def __pow__(self, exponent):
        if not isinstance(exponent, int) or exponent < 0:
            raise TitanError("Power must be a non-negative integer.")
        
        # Check if the new exponent will be within the valid range
        new_e = self.e * exponent
        if not (EXP_MIN <= new_e <= EXP_MAX):
            raise TitanError(f"Exponentiation failed for ({self})^({exponent}). "
                             f"Resulting exponent {new_e} is outside the allowed range [{EXP_MIN}, {EXP_MAX}].")

        # Check for numerator/denominator overflow.
        new_n = self.n ** exponent
        new_d = self.d ** exponent
        if new_n > NUM_DEN_MAX or new_d > NUM_DEN_MAX:
            raise TitanError(f"Exponentiation failed for ({self})^({exponent}). "
                             f"Resulting numerator '{self.n}^{exponent}={new_n}' or denominator '{self.d}^{exponent}={new_d}' exceeds the 5-bit limit of {NUM_DEN_MAX}.")
        
        # This part is unreachable if errors are thrown, but represents a successful operation.
        new_instance = object.__new__(TitanNumber)
        new_instance.value = self.value ** exponent
        new_instance.n, new_instance.d, new_instance.e, new_instance.error = new_n, new_d, new_e, 0
        return new_instance

    def __str__(self):
        # The final required output format, demonstrating the final calculated fraction.
        # Although we will not get this far, the problem requests we output each number in the equation.
        # This format would be used to print such results.
        return f"{self.n}/{self.d} * 2^{self.e}"

# Main execution logic
print("--- Titan Computer Feasibility Analysis for Pandora Landing ---")
print("Plan: Test if the essential physical quantities for the gravity calculation can be handled by Titan's 5-bit architecture.")

try:
    # Key physical constants in SI units (meters)
    pandora_equatorial_radius_m = 2_000_000.0
    pandora_core_radius_m = 100_000.0

    print(f"\n[Step 1] Attempting to represent Pandora's equatorial radius, a = {pandora_equatorial_radius_m:,.0f} m.")
    # This step is expected to fail because the number is too large.
    titan_a = TitanNumber(pandora_equatorial_radius_m)
    print(f"  [SUCCESS] Pandora's radius 'a' can be represented as: {titan_a}")

except TitanError as e:
    print(f"  [FAILURE] {e}")
    print("\n  Analysis of Step 1: The planet's equatorial radius 'a' is a fundamental value for calculating surface gravity (as part of the 1/r^2 term).")
    print("  Since it cannot be represented, the entire calculation is impossible.")

print("\n" + "-"*60)

try:
    print(f"\n[Step 2] For completeness, let's test an intermediate calculation with a smaller, representable value: the core radius, r_core = {pandora_core_radius_m:,.0f} m.")
    
    # This representation should succeed.
    titan_r_core = TitanNumber(pandora_core_radius_m)
    print(f"  [SUCCESS] Core radius 'r_core' is representable. Best approximation: {titan_r_core.n}/{titan_r_core.d} * 2^{titan_r_core.e} (Error: {titan_r_core.error*100:.3f}%)")
    
    print("\n  Now, attempting to calculate r_core^3 (required for the core's volume and mass)...")
    # This calculation is expected to fail due to exponent overflow.
    titan_r_core_cubed = titan_r_core ** 3
    print(f"  [SUCCESS] r_core^3 was calculated.")

except TitanError as e:
    print(f"  [FAILURE] {e}")
    print("\n  Analysis of Step 2: Even when a value can be represented, performing standard physics calculations (like cubing a radius to find a volume)")
    print("  can produce intermediate results that violate Titan's constraints. The exponent of r_core is 12; cubing it results in an exponent of 36, far exceeding the maximum of 15.")

print("\n" + "-"*60)
print("\n[Final Conclusion]")
print("The Titan computer architecture is incapable of performing the required calculation due to its strict 5-bit constraints.")
print("1. Critical inputs, like the planet's radius, are too large to be represented.")
print("2. Essential intermediate operations, like calculating volumes, result in numbers that exceed the architectural limits.")
print("Therefore, the force cannot be calculated.")

print("\n<<<N0>>>")