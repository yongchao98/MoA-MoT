# Titan Architecture Simulation for Pandora Escape Velocity

# This script simulates the step-by-step calculation on the Titan computer.
# It uses a "Register" class to hold expressions as lists of fractions.
# The core logic is in the 'titan_mul' function, which respects the 4-bit constraints.

MAX_INT = 15

class Register:
    """A class to simulate a Titan register holding an expression."""
    def __init__(self, numerator=0, denominator=1):
        # An expression is a list of terms (num, den).
        # A single fraction is an expression with one term.
        if numerator > MAX_INT or denominator > MAX_INT:
            raise ValueError(f"Constraint Violation: {numerator}/{denominator} exceeds 4-bit limit.")
        self.expression = [(numerator, denominator)]

    def set_expression(self, expr_list):
        self.expression = expr_list

    def __str__(self):
        return ' + '.join([f"{n}/{d}" for n, d in self.expression])

def titan_mul(reg1, reg2):
    """
    Simulates the MUL instruction. Multiplies two registers (reg1 * reg2).
    Returns a new expression list or raises a ValueError on violation.
    """
    new_expression = []
    # Distribute the multiplication across all terms in both expressions.
    for n1, d1 in reg1.expression:
        for n2, d2 in reg2.expression:
            # Check for direct multiplication violation
            if n1 * n2 > MAX_INT:
                # This is the point of failure. The Titan architecture cannot proceed
                # when a multiplication results in a numerator > 15 that cannot
                # be simplified by a denominator.
                print(f"--- COMPUTATION FAILED ---")
                print(f"Reason: During multiplication of the expression '{reg1}'")
                print(f"by the term '{n2}/{d2}', the operation '{n1}/{d1} * {n2}/{d2}' was required.")
                print(f"The resulting numerator '{n1} * {n2} = {n1*n2}' exceeds the 4-bit limit of {MAX_INT}.")
                print("This operation is not possible on the Titan architecture.")
                raise ValueError("Constraint Violation")
            
            # Simplification: Find common factors before multiplying
            # (e.g., a/b * c/a = c/b)
            # This is the most generous interpretation of the rules.
            import math
            common_factor1 = math.gcd(n1, d2)
            common_factor2 = math.gcd(n2, d1)
            
            final_num = (n1 // common_factor1) * (n2 // common_factor2)
            final_den = (d1 // common_factor2) * (d2 // common_factor1)
            
            if final_num > MAX_INT or final_den > MAX_INT:
                print(f"--- COMPUTATION FAILED ---")
                print(f"Reason: Even after simplification, multiplying '{n1}/{d1} * {n2}/{d2}'")
                print(f"results in '{final_num}/{final_den}', which exceeds the 4-bit limit of {MAX_INT}.")
                raise ValueError("Constraint Violation")

            new_expression.append((final_num, final_den))
    return new_expression


def run_pandora_calculation():
    """Main simulation logic."""
    print("Task: Calculate Pandora's escape velocity coefficient.")
    print("Formula: v_e^2 = (8/3) * pi * G * rho * R^2")
    print("We will calculate the coefficient part: C = (8/3) * pi * G_c * rho_c * R_c^2\n")

    try:
        # Step 1: Initialize registers with chosen approximations.
        # AX will be our accumulator register.
        # Let's start with the constant 8/3 from the formula.
        print("MOV AX, 8/3")
        AX = Register(8, 3)
        print(f"AX = {AX}\n")

        # Step 2: Multiply by Pi (pi approx 13/4 = 3 + 1/4)
        print("MUL AX, pi; where pi is approximated as 13/4 and expanded to (3/1 + 1/4)")
        pi_reg = Register()
        pi_reg.set_expression([(3, 1), (1, 4)]) # 3 + 1/4
        AX.set_expression(titan_mul(AX, pi_reg))
        print(f"AX = {AX}  (Represents 8/3 * (3+1/4) = 8 + 2/3)\n")

        # Step 3: Multiply by G_coefficient (G_c approx 2/3)
        # Here we attempt to calculate (8 + 2/3) * (2/3)
        print("MUL AX, G_c; where G_c is approximated as 2/3")
        G_reg = Register(2, 3)
        AX.set_expression(titan_mul(AX, G_reg))
        # This will fail because the first term multiplication is 8/1 * 2/3,
        # which results in a numerator of 16, violating the constraint.
        # The script will raise an exception and print the failure reason.
        
        print(f"AX = {AX}\n")
        print("--- CALCULATION SUCCEEDED ---")


    except ValueError:
        print("\nConclusion: The calculation cannot be completed on Titan.")
        print("The chain of multiplications required leads to numbers that exceed")
        print("the 4-bit representation, even when applying expansion and simplification rules.")
        # This confirms that the final answer should be N0.

run_pandora_calculation()
print("\nFinal Answer: Since the calculation is not possible, the answer is N0.")
# The following line is for the final answer format.
# print("<<<N0>>>")