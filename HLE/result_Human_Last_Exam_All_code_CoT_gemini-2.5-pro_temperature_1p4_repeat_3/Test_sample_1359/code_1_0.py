def solve_puzzle():
    """
    Solves the matrix puzzle by calculating the missing elements and their sum.
    """
    # M(3,1) is the starting triplet for the horizontal calculation
    m31 = [7, 2, 9]
    print(f"Starting with the given triplet M(3,1) = {m31}\n")

    # --- Calculate the first missing triplet: M(3,2) ---
    print("Step 1: Calculate the first missing triplet M(3,2) from M(3,1)")
    x, y, z = m31
    print(f"Current triplet [x y z] is [{x} {y} {z}]")
    
    # Check the horizontal transformation condition
    condition_sum = x + y
    print(f"Checking condition: x + y = {x} + {y} = {condition_sum}")

    if condition_sum > 10:
        # This block is not executed in this case
        next_x = (x * 3 - y) % 12
        next_y = (y * 2 + 4) % 12
        next_z = (z + x) % 12
    else:
        # x + y <= 10
        print(f"Since {condition_sum} <= 10, applying the second set of horizontal rules:")
        next_x = (x * 2 + y) % 12
        print(f"Next x = ({x} * 2 + {y}) mod 12 = {x*2+y} mod 12 = {next_x}")
        next_y = (y * 3 - 2) % 12
        print(f"Next y = ({y} * 3 - {y}) mod 12 = {y*3-2} mod 12 = {next_y}")
        next_z = (z * 2) % 12
        print(f"Next z = ({z} * 2) mod 12 = {z*2} mod 12 = {next_z}")

    m32 = [next_x, next_y, next_z]
    print(f"\nThe first missing triplet M(3,2) is: {m32}\n")

    # --- Calculate the second missing triplet: M(3,3) ---
    print("Step 2: Calculate the second missing triplet M(3,3) from M(3,2)")
    x, y, z = m32
    print(f"Current triplet [x y z] is [{x} {y} {z}]")
    
    # Check the horizontal transformation condition
    condition_sum = x + y
    print(f"Checking condition: x + y = {x} + {y} = {condition_sum}")

    if condition_sum > 10:
        # This block is not executed in this case
        next_x = (x * 3 - y) % 12
        next_y = (y * 2 + 4) % 12
        next_z = (z + x) % 12
    else:
        # x + y <= 10
        print(f"Since {condition_sum} <= 10, applying the second set of horizontal rules:")
        next_x = (x * 2 + y) % 12
        print(f"Next x = ({x} * 2 + {y}) mod 12 = {x*2+y} mod 12 = {next_x}")
        next_y = (y * 3 - 2) % 12
        print(f"Next y = ({y} * 3 - {y}) mod 12 = {y*3-2} mod 12 = {next_y}")
        next_z = (z * 2) % 12
        print(f"Next z = ({z} * 2) mod 12 = {z*2} mod 12 = {next_z}")
    
    m33 = [next_x, next_y, next_z]
    print(f"\nThe second missing triplet M(3,3) is: {m33}\n")
    
    # --- Calculate the sum of all missing elements ---
    print("Step 3: Calculate the sum of all missing elements")
    missing_elements = m32 + m33
    total_sum = sum(missing_elements)
    
    # Create the equation string as requested
    equation_str = " + ".join(map(str, missing_elements))
    
    print(f"The missing elements are: {', '.join(map(str, missing_elements))}")
    print(f"Final Sum = {equation_str} = {total_sum}")
    
    # Final answer in the required format
    # print(f"\n<<<{total_sum}>>>")

if __name__ == '__main__':
    solve_puzzle()