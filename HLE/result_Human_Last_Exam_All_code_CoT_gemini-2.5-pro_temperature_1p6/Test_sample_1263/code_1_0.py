import math

class TitanRegister:
    """
    A simplified simulator for a Titan register holding a single fraction.
    This demonstrates how the 4-bit constraint is violated.
    """
    def __init__(self):
        self.num = 0
        self.den = 1
        self.exp = 0
        print("Register initialized.")

    def check_constraint(self):
        """Checks if the current register value is valid."""
        if not (0 <= self.num <= 15 and 0 <= self.den <= 15 and self.den != 0):
            raise ValueError(f"FATAL ERROR: Constraint Violated. Value: {self.num}/{self.den}")
        return True

    def mov(self, val_str: str):
        """Simulates the MOV instruction."""
        print(f"\nExecuting: MOV AX, {val_str}")
        if 'e' in val_str:
            frac_part, exp_part = val_str.split('e')
            self.exp = int(exp_part)
        else:
            frac_part = val_str
            self.exp = 0
        
        if '/' in frac_part:
            num, den = map(int, frac_part.split('/'))
            self.num = num
            self.den = den
        else:
            self.num = int(frac_part)
            self.den = 1
        
        print(f" -> Register holds: {self.num}/{self.den} e{self.exp}")
        self.check_constraint()

    def mul(self, val_str: str):
        """Simulates the MUL instruction."""
        print(f"\nExecuting: MUL AX, {val_str}")
        
        # Parse operand
        if 'e' in val_str:
            frac_part, exp_part = val_str.split('e')
            exp2 = int(exp_part)
        else:
            frac_part = val_str
            exp2 = 0
        
        if '/' in frac_part:
            num2, den2 = map(int, frac_part.split('/'))
        else:
            num2 = int(frac_part)
            den2 = 1

        if not (0 <= num2 <= 15 and 0 <= den2 <= 15 and den2 != 0):
             raise ValueError(f"Invalid operand: {val_str}")

        # Perform multiplication
        new_num = self.num * num2
        new_den = self.den * den2
        new_exp = self.exp + exp2
        print(f" -> Intermediate product: {new_num}/{new_den} e{new_exp}")
        
        # Simplify the fraction
        common_divisor = math.gcd(new_num, new_den)
        self.num = new_num // common_divisor
        self.den = new_den // common_divisor
        self.exp = new_exp
        
        print(f" -> Simplified to: {self.num}/{self.den} e{self.exp}")
        
        # Check constraint after simplification
        self.check_constraint()

def demonstrate_impossibility():
    """
    Shows that calculating the necessary products for the escape velocity
    is not possible on Titan.
    """
    ax = TitanRegister()

    print("\n--- ATTEMPT 1: Multiply Pi and G first ---")
    try:
        ax.mov("13/4")  # Approximation for Pi
        ax.mul("13/2")  # Approximation for G. This will fail.
    except ValueError as e:
        print(f" -> {e}")

    print("\n--- ATTEMPT 2: Multiply other terms first ---")
    try:
        # These are mantissas for R^2, rho_s, and the factor 8
        ax.mov("4")
        ax.mul("3")
        ax.mul("8") # This will fail
    except ValueError as e:
        print(f" -> {e}")

    print("\nConclusion: The calculation fails because intermediate products systematically")
    print("exceed the 4-bit constraint for numerators (0-15), making it impossible.")
    print("Thus, the answer is N0.")
    print("\n<<<N0>>>")

demonstrate_impossibility()