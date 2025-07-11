import sys

def is_prime(n):
    """Checks if a number is prime. For this problem, we only care about primes up to 11."""
    return n in [2, 3, 5, 7, 11]

def calculate_missing_elements():
    """
    This function determines the missing elements based on the puzzle's rules and calculates their sum.
    """

    # Given triplet T(3,1)
    t31 = [7, 2, 9]

    # --- Step 1: Calculate the middle triplet T(3,2) from T(3,1) ---
    x1, y1, z1 = t31
    t32 = [0, 0, 0]

    # Apply Horizontal Transformations (Left -> Middle)
    if x1 + y1 > 10:
        # This condition is not met for T(3,1) as 7+2 = 9
        pass
    else: # x1 + y1 <= 10
        t32[0] = (x1 * 2 + y1) % 12
        t32[1] = (y1 * 3 - 2) % 12
        t32[2] = (z1 * 2) % 12

    # --- Step 2: Calculate the right triplet T(3,3) from T(3,2) ---
    x2, y2, z2 = t32
    t33 = [0, 0, 0]

    # Apply Cross-Dependency rule for x-value
    # This rule was deduced from the examples provided.
    # if middle z is prime: next x = (y + 2) % 12
    # if middle z is not prime: next x = (y - 1) % 12
    if is_prime(z2):
        t33[0] = (y2 + 2) % 12
    else:
        t33[0] = (y2 - 1) % 12
    
    # Apply Horizontal Transformations (Middle -> Right) for y and z values
    if x2 + y2 > 10:
        # This condition is not met for T(3,2) as 4+4=8
        pass
    else: # x2 + y2 <= 10
        t33[1] = (y2 * 3 - 2) % 12
        t33[2] = (z2 * 2) % 12

    # --- Step 3: Sum the missing elements and print the result ---
    missing_elements = t32 + t33
    total_sum = sum(missing_elements)
    
    # Format the final output equation string
    equation = " + ".join(map(str, missing_elements))
    
    # Print the equation and the final result
    print("The calculated missing triplet in the middle is: [{0}, {1}, {2}]".format(t32[0], t32[1], t32[2]))
    print("The calculated missing triplet on the right is: [{0}, {1}, {2}]".format(t33[0], t33[1], t33[2]))
    print("\nThe sum of the missing elements is calculated as follows:")
    print("{0} = {1}".format(equation, total_sum))
    
    # The final answer in the specified format
    # Redirecting to a variable to be appended later, to avoid printing it here.
    # This allows keeping the helpful prints above.
    final_answer_output = f"<<<{total_sum}>>>"
    
    # Temporarily redirect stdout to capture just the final answer
    original_stdout = sys.stdout
    sys.stdout = sys.stderr  # Redirect to stderr to hide from the main output
    print(final_answer_output)
    sys.stdout = original_stdout # Restore stdout


calculate_missing_elements()
