import sys

def compute_limit_expression(k):
    """
    Computes the value of the limit, which is 1 - 1/(2k).

    Args:
        k (int): An integer k >= 2.

    Returns:
        float: The value of the limit.
    """
    if not isinstance(k, int) or k < 2:
        raise ValueError("k must be an integer greater than or equal to 2.")
    
    numerator = 2 * k - 1
    denominator = 2 * k
    
    return float(numerator) / denominator

def main():
    """
    Main function to get user input for k and print the result.
    """
    print("The problem is to compute the limit of ln(f(m))/ln(m) as m -> infinity.")
    print("The result of the limit is given by the formula: (2*k - 1) / (2*k)")
    print("-" * 20)
    
    try:
        # Prompt the user to enter a value for k
        k_input = input("Please enter an integer value for k (k >= 2): ")
        k = int(k_input)
        
        # Calculate the limit
        limit_value = compute_limit_expression(k)
        
        # As per the instruction "output each number in the final equation!",
        # we identify the numbers in the formula (2*k - 1)/(2*k) as 2 and 1.
        print("\nThe numbers in the final equation (2*k - 1)/(2*k) are: 2, 1")
        
        # Print the final result for the given k
        print(f"For k = {k}, the value of the limit is: {limit_value}")

    except (ValueError, TypeError) as e:
        print(f"Invalid input: {e}. Please run the script again and enter a valid integer for k.", file=sys.stderr)
    except EOFError:
        print("\nNo input received. Please run the script again and provide input.", file=sys.stderr)


if __name__ == "__main__":
    main()
