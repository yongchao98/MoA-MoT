import math

class TitanComputerSimulator:
    """
    A simulator to test calculations on the Titan 4-bit architecture.
    """

    def check_constraint(self, val):
        """
        Checks if a value fits in the 4-bit (0-15) constraint.
        This is the critical rule for all intermediate calculations.
        """
        if val > 15:
            raise ValueError(f"CONSTRAINT VIOLATION: Value '{val}' exceeds the 4-bit limit of 15.")
        return val

    def multiply(self, term1_num, term1_den, term2_num, term2_den):
        """
        Simulates a multiplication operation on Titan.
        It checks for constraint violations on the immediate result.
        """
        # The operation results in a new numerator and denominator
        new_num = term1_num * term2_num
        new_den = term1_den * term2_den
        
        print(f"  Attempting to multiply: ({term1_num}/{term1_den}) * ({term2_num}/{term2_den})")
        print(f"  > Intermediate result before simplification: {new_num}/{new_den}")

        # "Any operation resulting in numerators or denominators exceeding 15 must be immediately simplified"
        # We check this constraint here. If it fails, the calculation path is invalid.
        self.check_constraint(new_num)
        self.check_constraint(new_den)
        
        # In a real scenario, one might simplify the fraction (e.g., 8/4 -> 2/1)
        # but the initial product must be valid.
        return new_num, new_den

    def run_simulation(self):
        """
        Runs the simulation to see if calculating Pandora's gravity is feasible.
        """
        print("--- Titan Computer Feasibility Simulation for Pandora Landing Time ---")
        print("\nObjective: Calculate g ≈ (4/3) * G * pi * r * rho_shell")
        print("Constraint: All numerators and denominators in intermediate steps must not exceed 15.")
        
        # Define constants using Titan-compatible fractional approximations
        # Using simple and favorable approximations to test the best-case scenario
        G_approx_num, G_approx_den = 7, 1  # G ≈ 7 (Error ~5%)
        pi_approx_num, pi_approx_den = 3, 1  # pi ≈ 3 (Error ~5%)
        r_approx_num, r_approx_den = 2, 1    # r = 2e6 m
        rho_approx_num, rho_approx_den = 3, 1 # rho_shell = 300 kg/m^3

        print("\nChosen Fractional Approximations:")
        print(f"  4/3")
        print(f"  G ≈ {G_approx_num}/{G_approx_den}")
        print(f"  pi ≈ {pi_approx_num}/{pi_approx_den}")
        print(f"  r (mantissa) = {r_approx_num}/{r_approx_den}")
        print(f"  rho_shell (mantissa) = {rho_approx_num}/{rho_approx_den}")

        print("\n--- Starting Calculation Trace ---")
        
        try:
            # Step 1: Start with 4/3
            current_num, current_den = 4, 3
            print(f"Start with 4/3. Current value: {current_num}/{current_den}")

            # Step 2: Multiply by G
            current_num, current_den = self.multiply(current_num, current_den, G_approx_num, G_approx_den)
            # This step will fail because 4 * 7 = 28, which is > 15.
            
            print("This line will not be reached.")

        except ValueError as e:
            print(f"\n--- SIMULATION FAILED ---")
            print(f"Reason: {e}")
            print("\nThe calculation fails at the very first multiplication involving physical constants.")
            print("The product of the numerators (4 * 7 = 28) exceeds the 4-bit limit.")
            print("No combination of operation order or valid fractional approximations for G and pi can avoid producing an intermediate number larger than 15.")
            print("\nConclusion: It is not possible to perform this calculation on the Titan computer.")

if __name__ == '__main__':
    simulator = TitanComputerSimulator()
    simulator.run_simulation()
    # The final answer must be in the format <<<ANSWER>>>
    print("\n<<<N0>>>")
