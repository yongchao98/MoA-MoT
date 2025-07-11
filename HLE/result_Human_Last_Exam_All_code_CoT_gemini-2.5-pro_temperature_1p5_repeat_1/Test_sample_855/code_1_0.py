import math

# Define the maximum value for a 5-bit unsigned integer
MAX_VAL = 31

class TitanComputer:
    """
    A class to simulate calculations on the Titan computer to determine
    the feasibility of the required landing force calculation.
    """

    def gcd(self, a, b):
        """Helper to compute the greatest common divisor for simplification."""
        return math.gcd(a, b)

    def multiply(self, frac1, frac2, name):
        """
        Simulates a multiplication operation on the Titan architecture.
        A calculation is successful only if the simplified result fits
        within the 5-bit constraints (<= 31).
        """
        n1, d1 = frac1
        n2, d2 = frac2

        print(f"Attempting to compute '{name}': ({n1}/{d1}) * ({n2}/{d2})")
        
        # Check if initial operands are valid for demonstration purposes
        if n1 > MAX_VAL or d1 > MAX_VAL or n2 > MAX_VAL or d2 > MAX_VAL:
             print(f"  - ERROR: An initial operand exceeds the 5-bit limit of {MAX_VAL}.")
             return None

        # Perform the multiplication of the fractions
        n_raw = n1 * n2
        d_raw = d1 * d2

        # Simplify the resulting fraction by dividing by the GCD
        common = self.gcd(n_raw, d_raw)
        n_final = n_raw // common
        d_final = d_raw // common

        print(f"  - Raw result before simplification: {n_raw}/{d_raw}")
        print(f"  - Simplified result after GCD: {n_final}/{d_final}")

        # Check if the final, simplified numbers adhere to the 5-bit constraint
        if n_final > MAX_VAL or d_final > MAX_VAL:
            print(f"  - FAILURE: Resulting numerator ({n_final}) or denominator ({d_final}) is > {MAX_VAL}.")
            return None
        
        print(f"  - SUCCESS: Result is within 5-bit constraints.")
        return (n_final, d_final)

    def run_feasibility_check(self):
        """
        Executes the simulation to check if the landing force can be calculated,
        step-by-step.
        """
        print("--- Titan Feasibility Study for Pioneer Landing Force ---")
        print("\nThe total required rocket force is F_rocket = F_gravity + F_deceleration.")
        print("We must check if both components can be calculated under 5-bit constraints.\n")

        # --- Check 1: Deceleration Force ---
        # The force for deceleration is F_decel = m*v^2/(2d) = 50 * 300^2 / (2*5000) = 450 N.
        # To calculate this, Titan must compute 50 * 9. With scientific notation, the
        # essential multiplication is (5/1) * (9/1).
        print("[1] Checking F_deceleration component (requiring 5 * 9)")
        op1_decel = (5, 1)
        op2_decel = (9, 1)
        self.multiply(op1_decel, op2_decel, name="Core of F_decel: 5 * 9")
        print()

        # --- Check 2: Gravitational Force ---
        # The force of gravity, F_g, depends on Pandora's mass, which involves products of
        # physical constants like G and pi. Let's test just the first part of this calculation.
        # g is proportional to (4/3) * pi * ...
        print("[2] Checking F_gravity component (requiring (4/3) * pi)")
        
        # We use standard, accurate fractional approximations for the constants.
        # pi ≈ 3.14159 -> 22/7
        four_thirds_frac = (4, 3)
        pi_frac = (22, 7)
        self.multiply(four_thirds_frac, pi_frac, name="Part of g calculation: (4/3) * pi")
        print()
        
        # --- Final Conclusion ---
        print("--- Final Result of Feasibility Study ---")
        print("The simulation shows that even the most fundamental calculations required")
        print("for this problem fail the 5-bit constraint. Multiplying common physical")
        print("constants and simple integers results in numerators that are larger than 31,")
        print("and these cannot be simplified back into the valid range.")
        print("\nBecause essential components of the final equation cannot be computed,")
        print("it is not possible to calculate the required force using Titan's architecture.")


# Execute the feasibility simulation
titan_simulator = TitanComputer()
titan_simulator.run_feasibility_check()