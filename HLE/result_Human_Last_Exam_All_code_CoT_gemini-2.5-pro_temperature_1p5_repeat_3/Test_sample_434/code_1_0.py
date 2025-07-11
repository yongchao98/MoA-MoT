import math

def find_blowup_range():
    """
    Calculates the range of y(0) for a given x(0) > 1 that leads to a blow-up
    in the solution of the given system of differential equations.
    """
    try:
        x0_str = input("Enter the value for the initial condition x(0) (must be > 1): ")
        x0 = float(x0_str)

        if x0 <= 1:
            print("Error: The value of x(0) must be greater than 1.")
            return

        # The blow-up condition is derived as y(0)^2 < 2*x(0) + 1 - 3*x(0)^(2/3).
        # We need to calculate the value of the right-hand side of this inequality.
        
        term1 = 2 * x0
        term2 = 1
        term3 = 3 * (x0**(2/3))
        
        rhs = term1 + term2 - term3
        
        # Analytically, for x0 > 1, rhs is always positive. We check for robustness.
        if rhs < 0:
            print("An unexpected calculation error occurred. The term inside the square root is negative.")
            return
            
        boundary = math.sqrt(rhs)

        print("\nFor the given system of ODEs, with x(0) =", x0)
        print("the solution blows up if the initial condition y(0) satisfies the inequality:")
        
        # Outputting the numbers in the final equation as requested.
        print(f"\ny(0)^2 < 2*x(0) + 1 - 3*x(0)^(2/3)")
        print(f"y(0)^2 < {term1} + {term2} - {term3:.6f}")
        print(f"y(0)^2 < {rhs:.6f}")
        print(f"|y(0)| < sqrt({rhs:.6f})")
        print(f"|y(0)| < {boundary:.6f}")

        print("\nThis means y(0) must be in the open interval:")
        print(f"y(0) ∈ (-{boundary:.6f}, {boundary:.6f})")

    except ValueError:
        print("\nInvalid input. Please enter a valid number for x(0).")
    except Exception as e:
        print(f"\nAn error occurred: {e}")

if __name__ == '__main__':
    find_blowup_range()