def calculate_toric_code_gsd():
    """
    Calculates the ground space degeneracy (GSD) of the toric code
    with n smooth holes and m rough holes.
    """
    try:
        n_str = input("Enter the number of smooth holes (n): ")
        n = int(n_str)
        if n < 0:
            print("Number of holes cannot be negative.")
            return

        m_str = input("Enter the number of rough holes (m): ")
        m = int(m_str)
        if m < 0:
            print("Number of holes cannot be negative.")
            return

        exponent = n + m
        degeneracy = 2 ** exponent

        # The final print statement is formatted to clearly show the equation
        print(f"The Ground Space Degeneracy (GSD) is calculated as 2^(n + m).")
        print(f"GSD = 2^({n} + {m}) = 2^{exponent} = {degeneracy}")

    except ValueError:
        print("Invalid input. Please enter integer values for n and m.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    calculate_toric_code_gsd()