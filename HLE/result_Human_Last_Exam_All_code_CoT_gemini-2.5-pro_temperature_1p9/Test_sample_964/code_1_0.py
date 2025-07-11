import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

def solve_chemistry_problem():
    """
    This function determines the final barium salt after a series of chemical reactions.
    """
    # Step 1: Barium chloride and silver nitrate react.
    # BaCl2(aq) + 2AgNO3(aq) -> Ba(NO3)2(aq) + 2AgCl(s)
    # The soluble barium salt formed is Barium Nitrate.

    # Step 2 & 3: Ammonia is added and then evaporated.
    # Ammonia reacts with AgCl but not with Ba(NO3)2. Removing ammonia
    # reverses the reaction with AgCl.
    # The barium salt remains unchanged throughout these steps.

    barium_salt_name = "Barium Nitrate"
    barium_salt_formula = "Ba(NO3)2"

    print(f"The barium salt present in the flask after all reactions is {barium_salt_name}.")
    print(f"Its chemical formula is: Ba(NO3)2")

    # Per the instructions, outputting each number from the final formula Ba(NO3)2
    # The number of nitrate groups is 2.
    # The number of oxygen atoms in each nitrate group is 3.
    print(f"The numbers in the chemical formula are 3 and 2.")

solve_chemistry_problem()

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the output to the console
print(output)

# Extract final answer from the captured output
answer = "Barium Nitrate"
print(f'<<<{answer}>>>', file=sys.stderr)