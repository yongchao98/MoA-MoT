import math
import argparse

def calculate_l_of_d(d):
    """
    Calculates the value of l(d) based on the analysis that the complex
    formula for f(d, p, o) simplifies to -sqrt(d).
    
    Args:
        d (int): The dimension of the space, must be an integer >= 2.
    """
    if not isinstance(d, int) or d < 2:
        print("Error: Input 'd' must be an integer greater than or equal to 2.")
        return

    result = -math.sqrt(d)
    
    # The prompt asks to "output each number in the final equation!".
    # This is interpreted as printing the full statement of the result,
    # for example: l(4) = -2.0.
    
    print(f"The final equation is l(d) = -sqrt(d)")
    print(f"For d = {d}, the numbers in the equation are:")
    print(f"l({d})")
    print("=")
    print(f"{result}")


if __name__ == '__main__':
    # We set up an argument parser to allow the user to specify the value of d.
    parser = argparse.ArgumentParser(
        description="Calculate the value of l(d) from the given problem.",
        epilog="The problem is analyzed to have a simple solution l(d) = -sqrt(d) hidden in a complex formula."
    )
    parser.add_argument('d', type=int, help="The dimension 'd', an integer >= 2.")
    
    args = parser.parse_args()
    
    calculate_l_of_d(args.d)
