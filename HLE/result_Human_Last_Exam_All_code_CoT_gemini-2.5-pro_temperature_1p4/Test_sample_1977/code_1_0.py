import math

def calculate_correlation_norm(n: int):
    """
    Calculates the 1-norm of the correlation matrix T for the state Jn with an even integer n.

    The formula used is: 3 + 2^(n+1) - (2^(n+2) + 4) / (1 + 3^n)

    Args:
        n: An even non-negative integer.
    """
    if not isinstance(n, int) or n < 0 or n % 2 != 0:
        print("Error: Please provide a non-negative even integer for n.")
        return

    # Calculate the components of the formula
    val_2_np1 = 2**(n + 1)
    val_2_np2 = 2**(n + 2)
    val_3_n = 3**n
    
    numerator = val_2_np2 + 4
    denominator = 1 + val_3_n
    
    # Calculate the final result
    # The result is guaranteed to be an integer, so we can safely convert it.
    result = 3 + val_2_np1 - numerator / denominator
    final_result = int(round(result))

    # Print the detailed breakdown of the calculation
    print(f"For n = {n}:")
    print(f"The formula for the 1-norm of the correlation matrix T is: 3 + 2^(n+1) - (2^(n+2) + 4) / (1 + 3^n)")
    print("\nBreaking down the calculation:")
    print(f"n + 1 = {n+1}")
    print(f"n + 2 = {n+2}")
    print(f"2^(n+1) = 2^{n+1} = {val_2_np1}")
    print(f"2^(n+2) = 2^{n+2} = {val_2_np2}")
    print(f"3^n = 3^{n} = {val_3_n}")
    
    print("\nSubstituting these values into the formula:")
    print(f"Result = 3 + {val_2_np1} - ({val_2_np2} + 4) / (1 + {val_3_n})")
    print(f"Result = {3 + val_2_np1} - {numerator} / {denominator}")
    print(f"Result = {3 + val_2_np1} - {numerator / denominator}")
    print(f"Result = {result}")

    print(f"\nThe final integer result is: {final_result}")


if __name__ == '__main__':
    # You can change the value of n here.
    # The problem asks for an even n. Let's use n=2 as an example.
    even_n = 2
    calculate_correlation_norm(even_n)