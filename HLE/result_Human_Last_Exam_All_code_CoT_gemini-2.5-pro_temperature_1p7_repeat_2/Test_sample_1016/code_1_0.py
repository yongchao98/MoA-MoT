import math

def solve_schwarz_iterations():
    """
    Calculates the number of iterations needed for the Schwarz method to converge
    for the 1D wave equation up to a given time T.
    """
    try:
        # Get user input for the parameters
        T_str = input("Enter the final time T: ")
        T = float(T_str)

        c_str = input("Enter the wave propagation speed c: ")
        c = float(c_str)

        M_str = input("Enter the overlap size M (M = b - a): ")
        M = float(M_str)

        if T < 0 or c <= 0 or M <= 0:
            print("Error: T must be non-negative, and c and M must be positive.")
            return

        # Calculate the number of iterations using the formula N = ceil(T * c / (2 * M))
        # The number of iterations must be an integer, so we use math.ceil.
        iterations = math.ceil((T * c) / (2 * M))

        # Print the final equation with the numbers plugged in
        print("\n--- Calculation ---")
        print(f"The number of iterations N is calculated as: ceil(T * c / (2 * M))")
        print(f"N = ceil({T} * {c} / (2 * {M}))")
        print(f"N = ceil({T * c} / {2 * M})")
        print(f"N = ceil({(T * c) / (2 * M)})")
        print(f"N = {iterations}")

    except ValueError:
        print("Error: Invalid input. Please enter valid numerical values.")
    except ZeroDivisionError:
        print("Error: Overlap size M cannot be zero.")

if __name__ == "__main__":
    solve_schwarz_iterations()