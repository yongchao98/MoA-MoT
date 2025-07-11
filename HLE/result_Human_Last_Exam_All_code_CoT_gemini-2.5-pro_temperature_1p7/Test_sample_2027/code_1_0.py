import math
import sys

def solve():
    """
    Calculates the exact value of l(d).
    The problem, after being corrected for likely typos, simplifies to the calculation of
    l(d) = 1 - max B(p, o)
    where B is the Busemann function. The maximum of the Busemann function B(p, o)
    for a fixed p is known to be the hyperbolic distance from the origin to p, which
    is 2 * artanh(||p||).
    Given p = (1/d, ..., 1/d), the Euclidean norm ||p|| is 1/sqrt(d).
    Thus, l(d) = 1 - 2 * artanh(1/sqrt(d)).
    """
    try:
        # Prompt the user to enter the value of d
        d_str = input("Enter the value for d (an integer >= 2): ")
        d = int(d_str)

        if d < 2:
            print("Error: d must be an integer greater than or equal to 2.", file=sys.stderr)
            return

        # Calculate the components of the final formula
        one_val = 1
        two_val = 2
        one_over_sqrt_d = 1 / math.sqrt(d)

        # Calculate l(d)
        l_d = one_val - two_val * math.atanh(one_over_sqrt_d)

        # Output the equation with the numbers used
        print(f"\nThe equation for l(d) is: l(d) = 1 - 2 * artanh(1/sqrt(d))")
        print(f"For d = {d}, we calculate:")
        # The instruction "output each number in the final equation" is interpreted here
        # by showing the values substituted into the formula.
        print(f"l({d}) = {one_val} - {two_val} * artanh({one_val}/sqrt({d}))")
        print(f"l({d}) = {one_val} - {two_val} * artanh({one_over_sqrt_d:.8f})")
        print(f"l({d}) = {l_d}")

    except ValueError:
        print("Invalid input. Please enter an integer.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)

# To run the solver
if __name__ == "__main__":
    solve()