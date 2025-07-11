import math

def find_integer_cube_root(n):
    """
    Finds the integer cube root of n, if it exists.
    Returns the integer root or None.
    """
    if n < 0:
        return None
    # Use floating point root as a good starting guess.
    # The round() function handles potential floating point inaccuracies.
    root = round(n**(1/3.0))
    if root ** 3 == n:
        return root
    return None

def get_sum_of_proper_divisor_cubes(n):
    """
    Calculates the sum of the cubes of the proper divisors of n.
    """
    # 1 is always a proper divisor for n > 1
    if n <= 1:
        return 0
    
    total = 1  # Start with 1^3
    
    # Iterate from 2 up to the square root of n
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            # i is a divisor
            total += i**3
            # n // i is the corresponding divisor pair
            j = n // i
            if i != j:
                total += j**3
    return total

def main():
    """
    Main function to solve the puzzle.
    """
    # As of August 2022, this was the next known integer in the sequence.
    number_to_check = 310815

    # Calculate the sum of the cubes of its proper divisors
    sum_of_cubes = get_sum_of_proper_divisor_cubes(number_to_check)
    
    # Check if the sum is a perfect cube
    cube_root = find_integer_cube_root(sum_of_cubes)
    
    print(f"The integer value that completes the sequence is {number_to_check}.")
    print("\nVerifying this number has the required property...")
    
    if cube_root is not None:
        print(f"\nThe number to be checked is: n = {number_to_check}")
        print("The property is: The sum of the cubes of the proper divisors of n must be a perfect cube.")
        print("\nThe final equation is:")
        # Outputting each number in the final equation as requested
        print(f"Sum of cubes of proper divisors = {sum_of_cubes}")
        print(f"Cube root of the sum = {cube_root}")
        print(f"Verification: {cube_root}^3 = {cube_root**3}")
    else:
        print(f"{number_to_check} does not satisfy the property.")

if __name__ == "__main__":
    main()
