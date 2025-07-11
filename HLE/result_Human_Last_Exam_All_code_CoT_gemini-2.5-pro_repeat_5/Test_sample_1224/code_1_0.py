import math

# This class simulates a number in the Titan computer's fractional system.
# It enforces the 4-bit constraint where numerators and denominators
# must be between 0 and 15.
class TitanFraction:
    MAX_VAL = 15

    def __init__(self, numerator, denominator=1):
        if not (0 <= numerator <= self.MAX_VAL and 0 < denominator <= self.MAX_VAL):
            raise ValueError(f"Initialization failed: Numerator/Denominator {numerator}/{denominator} out of 4-bit range (0-{self.MAX_VAL}).")
        self.num = numerator
        self.den = denominator
        self._simplify()

    def _simplify(self):
        if self.num == 0:
            self.den = 1
            return
        common = math.gcd(self.num, self.den)
        self.num //= common
        self.den //= common

    def __repr__(self):
        return f"{self.num}/{self.den}"

    def __mul__(self, other):
        # Simulate the multiplication operation on Titan
        print(f"Executing: ({self}) * ({other})")
        
        # Intermediate products before simplification
        res_num = self.num * other.num
        res_den = self.den * other.den
        
        # Titan rule: "Any operation resulting in numerators or denominators
        # exceeding 15 must be immediately simplified."
        # We check if simplification can bring the result back into range.
        common = math.gcd(res_num, res_den)
        final_num = res_num // common
        final_den = res_den // common

        if final_num > self.MAX_VAL or final_den > self.MAX_VAL:
            print(f"  - Intermediate result: {res_num}/{res_den}")
            print(f"  - Simplified result: {final_num}/{final_den}")
            raise ValueError(f"FATAL ERROR: Constraint violation. Resulting numerator '{final_num}' or denominator '{final_den}' exceeds {self.MAX_VAL}.")
        
        print(f"  - Result after simplification: {final_num}/{final_den}. Operation successful.")
        return TitanFraction(final_num, final_den)

def calculate_landing_time():
    """
    This function attempts to calculate the gravitational acceleration of Pandora
    using the Titan computer architecture simulation.
    """
    try:
        print("--- Titan Computer Calculation Simulation ---")
        print("Objective: Calculate Pandora's gravitational acceleration 'g'.")
        print("A simplified formula for g is proportional to: (4/3) * pi * r * rho")
        print("We will use integer coefficients for the calculation.\n")

        # Define constants as Titan Fractions
        # Volume of a sphere formula contains the factor 4/3
        v_factor = TitanFraction(4, 3)
        print(f"Factor from sphere volume formula: {v_factor}")

        # Planet's radius (r_planet) is 2000 km. We use the coefficient '2'.
        r_coeff = TitanFraction(2)
        print(f"Planet radius coefficient: {r_coeff}")
        
        # Shell density (rho_shell) is 0.3 metric tons/m^3. We use the coefficient '3'.
        rho_coeff = TitanFraction(3)
        print(f"Shell density coefficient: {rho_coeff}")
        
        # Pi is approximated as 3.
        pi_approx = TitanFraction(3)
        print(f"Pi approximation: {pi_approx}\n")

        # --- Begin Step-by-Step Calculation ---
        print("Step 1: Combine radius and density coefficients.")
        # This corresponds to the term (r_planet * rho_shell)
        term1 = r_coeff * rho_coeff # 2/1 * 3/1 = 6/1. Valid.

        print("\nStep 2: Multiply by the sphere volume factor (4/3).")
        # This corresponds to (4/3) * (r_planet * rho_shell)
        term2 = v_factor * term1 # 4/3 * 6/1 = 24/3, simplifies to 8/1. Valid.

        print("\nStep 3: Multiply by the approximation for Pi.")
        # This corresponds to (4/3) * pi * (r_planet * rho_shell)
        final_g_coeff = term2 * pi_approx # 8/1 * 3/1 = 24/1. Fails.

        # If we reached here, it means it was somehow possible.
        # This part of the code will not be reached.
        print("\nCalculation of g's coefficient was successful.")

    except ValueError as e:
        print(f"\n--- CALCULATION HALTED ---")
        print(e)
        print("\nThe Titan 4-bit architecture cannot handle the multiplication of intermediate values (e.g., 8 * 3 = 24)")
        print("even with the simplest approximations for the physical constants.")
        print("Therefore, calculating the landing time is not feasible.")
        print("\n<<<N0>>>")

if __name__ == '__main__':
    calculate_landing_time()