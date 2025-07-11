import math

class TitanComputer:
    """
    A simulator for the Titan 5-bit fractional computer.
    """
    MAX_VAL = 31 # 5-bit integer range (0-31)

    def find_best_fraction(self, value):
        """Finds the best fractional approximation for a value within 5-bit constraints."""
        best_frac = (0, 1)
        min_error = float('inf')

        for d in range(1, self.MAX_VAL + 1):
            # For a given denominator d, the best numerator is round(value * d)
            n = int(round(value * d))

            if 0 <= n <= self.MAX_VAL:
                error = abs(n / d - value)
                if error < min_error:
                    min_error = error
                    best_frac = (n, d)
        return best_frac

    def multiply(self, frac1, frac2):
        """
        Performs multiplication under Titan's rules.
        If direct multiplication exceeds 5-bit limits, it re-approximates.
        """
        n1, d1 = frac1
        n2, d2 = frac2

        # Check for direct multiplication feasibility
        if n1 * n2 <= self.MAX_VAL and d1 * d2 <= self.MAX_VAL:
            # Check for cross-cancellation, python's gcd is handy here
            common_nd = math.gcd(n1, d2)
            common_dn = math.gcd(n2, d1)
            
            new_n = (n1 // common_nd) * (n2 // common_dn)
            new_d = (d1 // common_dn) * (d2 // common_nd)

            if new_n <= self.MAX_VAL and new_d <= self.MAX_VAL:
                return (new_n, new_d), False

        # If direct multiplication is not possible, calculate true value and find best approximation
        true_value = (n1 / d1) * (n2 / d2)
        approximated_frac = self.find_best_fraction(true_value)
        return approximated_frac, True

def solve_coconut_problem():
    """
    Solves the Curious Monkey and Coconut problem using the TitanComputer simulator.
    """
    titan = TitanComputer()

    print("--- Titan Computer Analysis: Coconut Targeting System ---")
    print("\nStep 1: Deriving the governing equation for the force F.")
    print("The physics of a projectile under constant propulsion (F) and gravity (mg) leads to:")
    print("F = 2 * sqrt(2) * m * g")
    print("The mass 'm' is calculated from rock's radius (0.005 m) and density (900,000 kg/m^3):")
    print("m = (3/20) * pi kg")
    print("Substituting m, the final expression for the force F is:")
    print("F = (3/10) * sqrt(2) * pi * g\n")

    # --- Define high-precision constants for calculating the 'true' value ---
    g_true = 9.8
    pi_true = math.pi
    sqrt2_true = math.sqrt(2)
    
    # Calculate the ideal target force
    f_true = (3/10) * sqrt2_true * pi_true * g_true
    print(f"Step 2: Calculate the ideal target force using high-precision values.")
    print(f"F_ideal = (3/10) * {sqrt2_true:.4f} * {pi_true:.4f} * {g_true:.1f} = {f_true:.4f} N\n")

    print("Step 3: Perform the calculation using Titan's 5-bit fractional arithmetic.")
    print("We will approximate constants and re-approximate intermediate results if they exceed 5-bit limits.\n")

    # --- Titan Calculation ---
    # Initial term
    current_frac = (3, 10)
    print(f"Initial Term: 3/10")

    # Multiply by sqrt(2)
    f_sqrt2 = titan.find_best_fraction(sqrt2_true)
    print(f"Multiplying by sqrt(2), approximated as {f_sqrt2[0]}/{f_sqrt2[1]}...")
    current_frac, approx1 = titan.multiply(current_frac, f_sqrt2)
    if approx1:
        print(f"  Result exceeded 5-bit limit. Re-approximating to best fraction: {current_frac[0]}/{current_frac[1]}")
    
    # Multiply by pi
    f_pi = titan.find_best_fraction(pi_true)
    print(f"Multiplying by pi, approximated as {f_pi[0]}/{f_pi[1]}...")
    current_frac, approx2 = titan.multiply(current_frac, f_pi)
    if approx2:
        print(f"  Result exceeded 5-bit limit. Re-approximating to best fraction: {current_frac[0]}/{current_frac[1]}")

    # Multiply by g
    f_g = titan.find_best_fraction(g_true)
    print(f"Multiplying by g, approximated as {f_g[0]}/{f_g[1]}...")
    final_frac, approx3 = titan.multiply(current_frac, f_g)
    if approx3:
        print(f"  Result exceeded 5-bit limit. Re-approximating to best fraction: {final_frac[0]}/{final_frac[1]}")

    f_titan = final_frac[0] / final_frac[1]
    
    print("\n--- Calculation Summary ---")
    print("The final equation for the force F, calculated on Titan, is:")
    print(f"F = {final_frac[0]} / {final_frac[1]}")
    print(f"This corresponds to a decimal value of {f_titan:.4f} N.")

    # Calculate absolute error
    error = abs(f_titan - f_true)
    
    # Check if the force is within the allowed range [9.9m, 10.1m]
    # This corresponds to a force range of [0.99*F_true, 1.01*F_true]
    f_min = 0.99 * f_true
    f_max = 1.01 * f_true
    
    print(f"\nThe acceptable force range to hit the coconut is [{f_min:.2f}, {f_max:.2f}] N.")
    
    if f_min <= f_titan <= f_max:
        print(f"The calculated force {f_titan:.2f} N is within the acceptable range.")
        print("The system can hit the target.")
        answer = f"Y{error:.3f}"
    else:
        print(f"The calculated force {f_titan:.2f} N is outside the acceptable range.")
        print("The system cannot guarantee a hit.")
        answer = "N0"
        
    print("\n--- Final Answer ---")
    print(f"Absolute error |F_titan - F_ideal| = |{f_titan:.4f} - {f_true:.4f}| = {error:.4f}")
    print(f"Smallest absolute error 'e', rounded to 0.001, is {error:.3f}.")
    
    # The final required output format
    print(f"\n<<<{answer}>>>")

if __name__ == '__main__':
    solve_coconut_problem()
