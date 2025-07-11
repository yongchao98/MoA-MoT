import sys

# It's good practice to encapsulate the logic in functions.

def calculate_next_horizontal(triplet):
    """
    Applies the horizontal transformation rules to a given triplet.
    """
    x, y, z = triplet
    
    print(f"Applying horizontal transformation to [{x} {y} {z}]...")
    
    # Check the condition x + y
    if x + y > 10:
        print(f"Condition: {x} + {y} = {x+y} > 10")
        next_x = (x * 3 - y) % 12
        next_y = (y * 2 + 4) % 12
        next_z = (z + x) % 12
        print(f"Next x = ({x} * 3 - {y}) mod 12 = {next_x}")
        print(f"Next y = ({y} * 2 + 4) mod 12 = {next_y}")
        print(f"Next z = ({z} + {x}) mod 12 = {next_z}")
    else:
        print(f"Condition: {x} + {y} = {x+y} <= 10")
        next_x = (x * 2 + y) % 12
        next_y = (y * 3 - 2) % 12
        next_z = (z * 2) % 12
        print(f"Next x = ({x} * 2 + {y}) mod 12 = {next_x}")
        print(f"Next y = ({y} * 3 - 2) mod 12 = {next_y}")
        print(f"Next z = ({z} * 2) mod 12 = {next_z}")
        
    new_triplet = [next_x, next_y, next_z]
    print(f"Resulting triplet: {new_triplet}\n")
    return new_triplet

def solve_puzzle():
    """
    Solves the puzzle by calculating the missing elements and their sum.
    """
    # The starting triplet for the third row is given.
    m31 = [7, 2, 9]
    print(f"Starting with the first triplet of the third row: {m31}\n")

    # Calculate the second triplet in the third row (M(3, 2))
    m32 = calculate_next_horizontal(m31)

    # Calculate the third triplet in the third row (M(3, 3))
    m33 = calculate_next_horizontal(m32)

    # The missing elements are the numbers in m32 and m33.
    missing_elements = m32 + m33
    
    # Calculate the sum
    total_sum = sum(missing_elements)

    # Print the final sum equation as requested
    sum_equation = " + ".join(map(str, missing_elements))
    print("The missing elements are:", ", ".join(map(str, missing_elements)))
    print("Calculating the sum of the missing elements:")
    print(f"{sum_equation} = {total_sum}")
    
    # Final answer in the specified format
    # Redirecting the final answer to a specific stream is not standard.
    # We will print it clearly at the end.
    sys.stdout.write(f'<<<{total_sum}>>>')


# Run the solver
solve_puzzle()