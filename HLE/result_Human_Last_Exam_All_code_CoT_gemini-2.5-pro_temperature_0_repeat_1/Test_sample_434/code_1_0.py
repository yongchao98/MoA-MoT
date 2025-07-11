import math

def find_blowup_condition():
    """
    This function calculates the critical value of y(0) for a given x(0) > 1
    that determines whether the solution to a system of ODEs blows up.

    The system is:
    x'(t) = -3*x(t)*y(t)
    y'(t) = -y(t)**2 - x(t) + 1
    """
    print("This script calculates the range of y(0) for which the solution to the system blows up, given x(0) > 1.")
    print("The blow-up occurs if y(0) is less than a critical value determined by x(0).")
    print("The critical value is given by the expression: sqrt(2*x(0) + 1 - 3*x(0)^(2/3))")
    print("-" * 70)

    try:
        x0_str = input("Please enter the initial value x(0) (must be a number greater than 1): ")
        x0 = float(x0_str)

        if x0 <= 1:
            print("Error: The value of x(0) must be strictly greater than 1.")
            return

        # The expression under the square root is 2*x0 + 1 - 3*x0^(2/3)
        # Let's calculate each part as requested by the prompt.
        c1 = 2
        c2 = 1
        c3 = -3
        p = 2/3

        term1 = c1 * x0
        term2 = c2
        term3 = c3 * (x0**p)

        val_under_sqrt = term1 + term2 + term3

        # For x0 > 1, val_under_sqrt is always positive.
        critical_y0 = math.sqrt(val_under_sqrt)

        print(f"\nFor the initial condition x(0) = {x0}, the solution blows up if:")
        print(f"y(0) < sqrt({c1}*x(0) + {c2} + {c3}*x(0)^({p:.2f}))")
        
        print("\nSubstituting the value of x(0):")
        print(f"y(0) < sqrt({c1}*{x0} + {c2} + {c3}*({x0})^({p:.2f}))")
        print(f"y(0) < sqrt({term1} + {term2} {term3:+.4f})")
        print(f"y(0) < sqrt({val_under_sqrt:.4f})")
        print(f"y(0) < {critical_y0:.4f}")

        print("\nConclusion:")
        print(f"The solution blows up for any initial value y(0) in the interval (-infinity, {critical_y0:.4f}).")

    except ValueError:
        print("Invalid input. Please enter a valid number for x(0).")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    find_blowup_condition()