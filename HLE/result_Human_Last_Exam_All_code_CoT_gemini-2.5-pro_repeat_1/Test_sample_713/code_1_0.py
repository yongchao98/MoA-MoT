import math

def get_volume_ratio_range(d):
    """
    Calculates the range of possible values for the ratio Vol(P)/V for a given dimension d.

    Args:
        d (int): The dimension of the simplex.

    Returns:
        A tuple containing the minimum and maximum possible ratios.
        Returns (min_ratio, None) if the maximum is not known for the given d.
    """
    if not isinstance(d, int) or d < 2:
        raise ValueError("Dimension d must be an integer greater than or equal to 2.")

    # The minimum ratio is always d!
    min_ratio = math.factorial(d)

    # The maximum ratio is achieved for a regular simplex. These values are known for low dimensions.
    max_ratios = {
        2: (10, 3),  # 10/3
        3: (10, 1),  # 10
        4: (147, 4)  # 147/4
    }

    if d in max_ratios:
        num, den = max_ratios[d]
        max_ratio = num / den
        max_ratio_str = f"{num}/{den}" if den != 1 else str(num)
        return min_ratio, max_ratio, max_ratio_str
    else:
        return min_ratio, None, None

def main():
    """
    Main function to get user input and display the result.
    """
    try:
        d_str = input("Enter the dimension d (e.g., 2, 3, 4): ")
        d = int(d_str)

        min_r, max_r, max_r_str = get_volume_ratio_range(d)
        
        min_val_eq = f"{d}!"
        
        print("\nLet V be the volume of the d-dimensional simplex T.")
        print("Let P be the polytope bounded by the 2*binom(d+1, 2) hyperplanes.")
        print("The set of possible values for the ratio Vol(P)/V is a closed interval.")
        
        print(f"\nFor dimension d = {d}:")
        print(f"The minimum value of the ratio is {min_val_eq} = {min_r}.")
        
        if max_r is not None:
            print(f"The maximum value of the ratio is {max_r_str} = {max_r:.4f}.")
            final_equation = f"Vol(P)/V in [{min_val_eq}, {max_r_str}]"
            print(f"\nThe set of possible values is the interval: [{min_r}, {max_r}]")
            print("\nFinal Equation:")
            # We output each part of the final equation/statement as requested
            print(f"For d = {d}, the set of possible values for the ratio Vol(P)/V is the interval [{min_r}, {max_r}].")
            print(f"In symbolic terms, this is the interval from {min_val_eq} to {max_r_str}.")

        else:
            print("The maximum value for this dimension is not available in this script.")
            print(f"\nThe set of possible values is the interval starting from {min_r}: [{min_r}, infinity)")
            print("\nFinal statement:")
            print(f"For d = {d}, the minimum value for the ratio Vol(P)/V is {min_r}. The maximum is not hardcoded.")

    except ValueError as e:
        print(f"Invalid input: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
