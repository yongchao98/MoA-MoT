import math

def calculate_area_ratio(n):
    """
    Calculates the ratio of the area of a constructed n-sided regular polygon
    to the starting 2n-sided regular polygon.

    The general formula for this ratio is: (1 + sec(pi/n)) / 2
    """
    if not isinstance(n, int) or n <= 2:
        print(f"Error: n must be an integer greater than 2.")
        return

    # Calculate the angle pi/n in radians
    angle_in_radians = math.pi / n

    # Calculate sec(pi/n). math.sec doesn't exist, so we use 1/cos().
    sec_value = 1 / math.cos(angle_in_radians)

    # Calculate the final ratio using the formula
    ratio = (1 + sec_value) / 2

    # As requested, output the numbers in the final equation
    print(f"For n = {n} (starting with a 2n = {2*n} sided polygon):")
    print("The formula is: (1 + sec(pi/n)) / 2")
    print(f"The calculation is: (1 + sec(pi/{n})) / 2")
    print(f"= (1 + {sec_value:.4f}) / 2")
    print(f"= {1 + sec_value:.4f} / 2")
    print(f"= {ratio:.4f}")
    print(f"The area of the {n}-sided polygon is {ratio:.4f} times the area of the {2*n}-sided polygon.")

# Main execution: Solve for the n=3 case from the problem description
if __name__ == "__main__":
    # Case from the prompt: 2n=6 sided hexagon, n=3 sided triangle
    n_value = 3
    calculate_area_ratio(n_value)
