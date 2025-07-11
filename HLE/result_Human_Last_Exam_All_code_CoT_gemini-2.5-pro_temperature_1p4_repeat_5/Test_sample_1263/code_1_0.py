import math

# A helper class to represent a single fraction with an exponent (A/B * 10^E)
# This corresponds to a single term in a Titan register.
class Term:
    def __init__(self, num, den, exp=0):
        if not (isinstance(num, int) and isinstance(den, int) and 0 <= num <= 15 and 1 <= den <= 15):
            raise ValueError(f"FATAL: Invalid 4-bit integer in constructor: {num}/{den}")
        self.num = num
        self.den = den
        self.exp = exp

    def __repr__(self):
        # This function helps in printing the term as per the problem's notation.
        # e.g., 12/1e14
        return f"{self.num}/{self.den}e{self.exp}"

# A helper function to print the equation in the register.
def print_equation(register_name, expression):
    # This function prints each number in the equation stored in the register.
    if not expression.terms:
        print(f"Register {register_name} is empty (0)")
        return
        
    equation_parts = []
    for term in expression.terms:
        equation_parts.append(f"{term.num}/{term.den}e{term.exp}")
    
    print(f"Register {register_name} = " + " + ".join(equation_parts))


# A class to simulate a Titan register, which holds an expression (a sum of terms).
class TitanRegister:
    def __init__(self, name="AX"):
        self.name = name
        self.terms = []
        print(f"Initialized register {self.name}")

    def mov(self, term):
        """Simulates the MOV instruction."""
        print(f"MOV {self.name}, {term}")
        self.terms = [term]
        print_equation(self.name, self)

    def add(self, term):
        """Simulates the ADD instruction."""
        print(f"ADD {self.name}, {term}")
        self.terms.append(term)
        if len(self.terms) > 10:
            raise ValueError(f"Term limit exceeded: {len(self.terms)} terms.")
        print_equation(self.name, self)

    def mul(self, factor):
        """Simulates the MUL instruction with immediate simplification."""
        print(f"MUL {self.name}, {factor}")
        
        new_terms_list = []
        for term in self.terms:
            num = term.num * factor.num
            den = term.den * factor.den
            exp = term.exp + factor.exp
            
            # Simplify fraction by dividing by greatest common divisor
            common_divisor = math.gcd(num, den)
            num //= common_divisor
            den //= common_divisor
            
            # "Immediate simplification" for numbers > 15
            if num > 15 or den > 15:
                # For this problem, we only encounter numerators > 15
                if den > 15:
                    raise ValueError(f"Denominator {den} exceeded 15.")

                # Decompose the integer part
                integer_part = num // den
                rem_num = num % den
                
                while integer_part > 0:
                    val = min(integer_part, 15)
                    new_terms_list.append(Term(val, 1, exp))
                    integer_part -= val
                
                # Add the fractional part if it exists
                if rem_num > 0:
                    common_rem = math.gcd(rem_num, den)
                    new_terms_list.append(Term(rem_num // common_rem, den // common_rem, exp))
            else:
                new_terms_list.append(Term(num, den, exp))

        self.terms = new_terms_list
        print_equation(self.name, self)
        if len(self.terms) > 10:
            raise ValueError(f"Term limit exceeded after MUL: {len(self.terms)} terms.")

def calculate_escape_velocity():
    """
    Main function to simulate the calculation on Titan.
    """
    print("--- Starting Pandora Escape Velocity Calculation on Titan Simulator ---\n")
    try:
        # Initialize register AX
        ax = TitanRegister("AX")

        # Step 1: Calculate the main components of the mass term
        # Term A = rho_s * R_total^2 = (3e2) * (2e6)^2 = 12e14
        term_A = Term(12, 1, 14)
        # Term B = (rho_c - rho_s) * r_c^3 / R_total = (9e2) * (1e5)^3 / (2e6) = 4.5e11 = 9/2e11
        term_B = Term(9, 2, 11)
        
        ax.mov(term_A)
        ax.add(term_B)
        
        # Step 2: Sequentially multiply by constants for v_e^2
        # Factor C = 8/3 * G * pi
        # We use approximations: G ≈ 20/3e-11 = 10/1 * 2/3 * e-11, pi ≈ 22/7
        
        # Multiply by 8/3
        ax.mul(Term(8, 3, 0))

        # Multiply by G_approx part 1 (10/1)
        # This is where the term count is expected to explode and fail.
        ax.mul(Term(10, 1, 0))

        # The following steps would be executed if the previous one hadn't failed.
        # ax.mul(Term(2, 3, -11)) # G_approx part 2
        # ax.mul(Term(22, 7, 0))  # Pi_approx
        
        print("\n--- Calculation finished successfully (hypothetical) ---")

    except ValueError as e:
        print(f"\n--- CALCULATION FAILED ---")
        print(f"Error: {e}")
        print("The calculation is not feasible on Titan's architecture.")
        print("\nFinal Answer: N0")

if __name__ == '__main__':
    calculate_escape_velocity()