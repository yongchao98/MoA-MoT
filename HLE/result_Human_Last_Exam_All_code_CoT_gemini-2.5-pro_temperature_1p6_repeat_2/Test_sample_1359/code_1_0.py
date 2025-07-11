def solve_puzzle():
    """
    This function solves the puzzle by calculating the missing matrix elements
    and their sum based on the provided transformation rules.
    """

    # Helper function for modulo 12 arithmetic
    def mod12(n):
        return n % 12

    # Horizontal Transformation Rules
    def horizontal_transform(triplet):
        x, y, z = triplet
        if x + y > 10:
            # Rule H1
            next_x = mod12(x * 3 - y)
            next_y = mod12(y * 2 + 4)
            next_z = mod12(z + x)
        else:
            # Rule H2
            next_x = mod12(x * 2 + y)
            next_y = mod12(y * 3 - 2)
            next_z = mod12(z * 2)
        return [next_x, next_y, next_z]

    # Given known matrix triplet M[3][1]
    m31 = [7, 2, 9]
    print(f"Given the starting triplet in the last row: {m31}\n")

    # Step 1: Calculate the first missing triplet, M[3][2]
    # This is a Left -> Middle transformation, so we apply the horizontal rule to M[3][1].
    # For M[3][1] = [7, 2, 9], x+y = 7+2 = 9, which is not > 10. Rule H2 applies.
    print("Calculating the first missing triplet [? ? ?] from [7 2 9]:")
    x, y, z = m31
    # next_x = (x * 2 + y) mod 12 = (7 * 2 + 2) mod 12 = 16 mod 12 = 4
    # next_y = (y * 3 - 2) mod 12 = (2 * 3 - 2) mod 12 = 4 mod 12 = 4
    # next_z = (z * 2) mod 12 = (9 * 2) mod 12 = 18 mod 12 = 6
    m32 = horizontal_transform(m31)
    print(f"Applying rule x+y <= 10: Next x = ({x} * 2 + {y}) mod 12 = {m32[0]}")
    print(f"Applying rule x+y <= 10: Next y = ({y} * 3 - 2) mod 12 = {m32[1]}")
    print(f"Applying rule x+y <= 10: Next z = ({z} * 2) mod 12 = {m32[2]}")
    print(f"The first missing triplet is: {m32}\n")


    # Step 2: Calculate the second missing triplet, M[3][3]
    # This is a Middle -> Right transformation, apply horizontal rule to M[3][2].
    # For M[3][2] = [4, 4, 6], x+y = 4+4 = 8, which is not > 10. Rule H2 applies again.
    print("Calculating the second missing triplet [? ? ?] from [4 4 6]:")
    x, y, z = m32
    # next_x = (x * 2 + y) mod 12 = (4 * 2 + 4) mod 12 = 12 mod 12 = 0
    # next_y = (y * 3 - 2) mod 12 = (4 * 3 - 2) mod 12 = 10 mod 12 = 10
    # next_z = (z * 2) mod 12 = (6 * 2) mod 12 = 12 mod 12 = 0
    m33 = horizontal_transform(m32)
    print(f"Applying rule x+y <= 10: Next x = ({x} * 2 + {y}) mod 12 = {m33[0]}")
    print(f"Applying rule x+y <= 10: Next y = ({y} * 3 - 2) mod 12 = {m33[1]}")
    print(f"Applying rule x+y <= 10: Next z = ({z} * 2) mod 12 = {m33[2]}")
    print(f"The second missing triplet is: {m33}\n")


    # Step 3: Sum the missing elements
    missing_elements = m32 + m33
    total_sum = sum(missing_elements)
    
    print(f"The missing elements are: {', '.join(map(str, missing_elements))}")
    equation = " + ".join(map(str, missing_elements))
    print(f"The sum of the missing elements is: {equation} = {total_sum}")

solve_puzzle()
<<<24>>>