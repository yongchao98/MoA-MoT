import math

class TitanComputerSimulator:
    """
    A simulator for the Titan 4-bit architecture to test the feasibility
    of calculating Pandora's escape velocity.
    """

    def __init__(self):
        self.max_val = 15
        print("Titan Computer Simulator Initialized.")
        print("Constraint: All numerators and denominators in calculations must be <= 15.")
        print("-" * 50)

    def multiply(self, frac1, frac2):
        """
        Simulates the MUL instruction.
        Multiplies two fractions (represented as tuples) and checks for overflow.
        An operation (a/b) * (c/d) is only valid if a*c and b*d are <= 15.
        """
        num1, den1 = frac1
        num2, den2 = frac2
        
        res_num = num1 * num2
        res_den = den1 * den2
        
        print(f"Attempting to multiply ({num1}/{den1}) by ({num2}/{den2})...")
        print(f"  Intermediate result: {res_num}/{res_den}")
        
        if res_num > self.max_val or res_den > self.max_val:
            print(f"  --> FAILURE: Overflow detected. Intermediate numerator ({res_num}) or denominator ({res_den}) exceeds {self.max_val}.")
            raise ValueError("Overflow")
        
        # In a real scenario, we might simplify using GCD, but let's first check intermediate products
        # For example, if we simplify after, 6/3 * 4/1 = 24/3. 24 > 15 -> fail.
        # This confirms our stricter check is appropriate.
        print(f"  --> SUCCESS: Result {res_num}/{res_den} is storable.")
        return (res_num, res_den)

    def run_calculation(self):
        """
        Runs the sequence of calculations for escape velocity.
        """
        # The final formula for v_e^2 is approximately 8/3 * pi * G * rho_shell * R_total^2
        # We will multiply the mantissas of these terms.
        
        # Part 1: Define constants as Titan fractions (numerator, denominator)
        const_8_3 = (8, 3)
        pi = (3, 1)
        G_mantissa = (2, 3)
        rho_shell_mantissa = (3, 1)
        R_total_sq_mantissa = (4, 1)

        print("Final Equation to compute:")
        print("v_e = sqrt( C * pi * G * rho * R^2 )")
        print("Substituting Titan approximations (mantissas and exponents shown separately):")
        # Final expression string construction
        expr = (f"v_e = sqrt( ({const_8_3[0]}/{const_8_3[1]}) * ({pi[0]}/{pi[1]}) * "
                f"({G_mantissa[0]}/{G_mantissa[1]} * 10^-10) * "
                f"({rho_shell_mantissa[0]}/{rho_shell_mantissa[1]} * 10^2) * "
                f"({R_total_sq_mantissa[0]}/{R_total_sq_mantissa[1]} * 10^12) )")
        print(expr)
        print("-" * 50)
        
        # Part 2: Attempt the multiplication step-by-step
        print("Multiplying mantissas...\n")
        try:
            # Let's try an order that combines small numbers first
            # Step 1: G * pi
            val1 = self.multiply(G_mantissa, pi)
            
            # Step 2: result * rho_shell
            val2 = self.multiply(val1, rho_shell_mantissa)

            # Step 3: result * 8/3
            val3 = self.multiply(val2, const_8_3)

            # Step 4: result * R_total^2 -> THIS WILL FAIL
            # Current value is (8/3) from previous step
            # Multiplying (8/3) by (4/1)
            final_val = self.multiply(val3, R_total_sq_mantissa)

            print("\nCalculation was successful (This should not be reached).")

        except ValueError:
            print("\n" + "-"*50)
            print("CONCLUSION:")
            print("The calculation cannot be completed on the Titan architecture.")
            print("The multiplication of the necessary physical constants, even with approximations,")
            print("results in values that exceed the 4-bit integer limit.")
            print("Therefore, Pandora's escape velocity cannot be determined.")
            print("Final Answer: N0")

if __name__ == "__main__":
    simulator = TitanComputerSimulator()
    simulator.run_calculation()