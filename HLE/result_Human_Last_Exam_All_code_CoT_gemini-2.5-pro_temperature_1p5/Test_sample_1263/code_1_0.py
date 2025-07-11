import math

class TitanConstraintError(Exception):
    """Custom exception for Titan architecture violations."""
    def __init__(self, message):
        super().__init__(message)

class TitanNumber:
    """
    Represents a number in Titan's memory. It is stored as a list of terms,
    simulating an expression in a register. Each term is of the form (num/den)*10^exp.
    """
    def __init__(self, terms):
        self.terms = []
        for num, den, exp in terms:
            if not (0 <= num <= 15 and 1 <= den <= 15):
                raise TitanConstraintError(f"Invalid term: Numerator or denominator in {num}/{den} is outside the 4-bit range [0, 15].")
            self.terms.append({'num': num, 'den': den, 'exp': exp})

    @staticmethod
    def gcd(a, b):
        """Helper to calculate the greatest common divisor."""
        return math.gcd(a, b)

    def __mul__(self, other):
        """
        Overloads the multiplication operator to simulate the MUL instruction.
        It checks if the resulting fraction's components exceed the 4-bit limit.
        """
        result_terms = []
        # Multiplication distributes over addition: (t1+t2)*(t3+t4) = t1t3+t1t4+t2t3+t2t4
        for t1 in self.terms:
            for t2 in other.terms:
                # Perform the multiplication of the two fractional terms
                raw_num = t1['num'] * t2['num']
                raw_den = t1['den'] * t2['den']

                # Immediately simplify the resulting fraction using GCD
                common_divisor = self.gcd(raw_num, raw_den)
                final_num = raw_num // common_divisor
                final_den = raw_den // common_divisor

                # Check if the simplified result violates the 4-bit constraint
                if final_num > 15 or final_den > 15:
                    raise TitanConstraintError(
                        f"Multiplication failed for ({t1['num']}/{t1['den']}) * ({t2['num']}/{t2['den']}). "
                        f"The simplified result {final_num}/{final_den} has components exceeding 15."
                    )
                
                new_exp = t1['exp'] + t2['exp']
                result_terms.append({'num': final_num, 'den': final_den, 'exp': new_exp})
        
        if len(result_terms) > 10:
            raise TitanConstraintError("Expression exceeds 10 terms in register after multiplication.")
        return TitanNumber([(t['num'], t['den'], t['exp']) for t in result_terms])

def run_titan_calculation():
    """
    Simulates the calculation of Pandora's escape velocity on the Titan architecture.
    """
    try:
        # We need to compute v_e^2 ≈ (8/3) * pi * G * d_shell * R^2
        # Let's define the constants and variables as TitanNumbers.
        # We choose approximations with small integer components to maximize success chance.
        
        # G ≈ 6.674e-11. Let's use 13/2 = 6.5.
        G = TitanNumber([(13, 2, -11)])
        # pi ≈ 3.14159. Let's use a simple 3/1.
        PI = TitanNumber([(3, 1, 0)])
        # Planet's shell density d_shell = 300 kg/m^3 = 3e2
        D_SHELL = TitanNumber([(3, 1, 2)])
        # Planet's radius R = 2e6 m, so R^2 = 4e12
        R_SQUARED = TitanNumber([(4, 1, 12)])
        # The constant 8/3 from the formula
        EIGHT_THIRDS = TitanNumber([(8, 3, 0)])

        print("The target expression to calculate is v_e^2, which is approximately the product of:")
        print(f"({EIGHT_THIRDS.terms[0]['num']}/{EIGHT_THIRDS.terms[0]['den']}) * ({PI.terms[0]['num']}/{PI.terms[0]['den']}) * ({G.terms[0]['num']}/{G.terms[0]['den']}e{G.terms[0]['exp']}) * ({D_SHELL.terms[0]['num']}/{D_SHELL.terms[0]['den']}e{D_SHELL.terms[0]['exp']}) * ({R_SQUARED.terms[0]['num']}/{R_SQUARED.terms[0]['den']}e{R_SQUARED.terms[0]['exp']})")
        print("\n--- Starting Titan Calculation Simulation ---")
        
        # We will try a specific multiplication order to see if it's possible.
        # Let's start with simple numbers.
        
        # 1. MOV AX, 8/3; MUL AX, pi
        print(f"1. Calculating (8/3) * pi = ({EIGHT_THIRDS.terms[0]['num']}/{EIGHT_THIRDS.terms[0]['den']}) * ({PI.terms[0]['num']}/{PI.terms[0]['den']})")
        AX = EIGHT_THIRDS * PI
        print(f"   Intermediate Result: {AX.terms[0]['num']}/{AX.terms[0]['den']}")

        # 2. MUL AX, d_shell
        print(f"\n2. Multiplying by d_shell = ({AX.terms[0]['num']}/{AX.terms[0]['den']}) * ({D_SHELL.terms[0]['num']}/{D_SHELL.terms[0]['den']})")
        AX = AX * D_SHELL
        print(f"   Intermediate Result: {AX.terms[0]['num']}/{AX.terms[0]['den']}")
        
    except TitanConstraintError as e:
        print(f"\nERROR: A constraint of the Titan architecture was violated.")
        print(e)
        print("\nConclusion: The calculation is not feasible as it requires multiplying numbers that result in components larger than 15.")
        print("<<<N0>>>")

# Run the simulation
run_titan_calculation()