import math

def find_blowup_condition():
    """
    This function calculates the range of initial values for y(0)
    that lead to a blow-up solution for the given system of ODEs,
    assuming x(0) > 1.
    """
    try:
        x0_str = input("Enter the initial condition x(0) > 1: ")
        x0 = float(x0_str)

        if x0 <= 1:
            print("Error: The initial condition x(0) must be greater than 1.")
            return

        # The separatrix y_c(x) = sqrt(1 + 2*x - 3*x^(2/3)) defines the boundary.
        # A solution blows up if y(0) < y_c(x(0)).
        
        # Calculate each term of the equation
        term1 = 1
        term2 = 2 * x0
        x0_pow_2_3 = x0**(2/3)
        term3 = -3 * x0_pow_2_3
        
        # Calculate the value under the square root
        yc_squared = term1 + term2 + term3
        
        print(f"\nFor x(0) = {x0}:")
        print("The critical value for y(0) is derived from the separatrix equation:")
        print("y_c = sqrt(1 + 2*x(0) - 3*x(0)^(2/3))")
        
        print("\nCalculating the terms:")
        print(f"y_c^2 = 1 + 2*({x0}) - 3*({x0})^(2/3)")
        print(f"y_c^2 = {term1} + {term2} - 3*({x0_pow_2_3})")
        print(f"y_c^2 = {term1 + term2} {term3}")
        
        if yc_squared < 0:
            # This shouldn't happen for x0 > 1 as shown in the derivation
            print(f"y_c^2 = {yc_squared}")
            print("The value under the square root is negative, which indicates an issue.")
        else:
            yc = math.sqrt(yc_squared)
            print(f"y_c^2 = {yc_squared}")
            print(f"y_c = {yc}")
            print("\nTherefore, the solution blows up if y(0) is less than this critical value.")
            print(f"Blow-up condition: y(0) < {yc}")

    except ValueError:
        print("Invalid input. Please enter a number.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    find_blowup_condition()