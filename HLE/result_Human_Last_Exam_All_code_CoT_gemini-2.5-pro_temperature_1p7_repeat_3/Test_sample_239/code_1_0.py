# This script calculates the minimum number of flagstones required to pave a rectangular area.
# It uses an efficient, no-branching formula to find the solution.

def solve_theatre_square():
    """
    Reads the dimensions of the square (n, m) and the flagstone (a),
    and prints the number of flagstones needed in an equation format.
    """
    try:
        # 1. Read the input values for n, m, and a from the user.
        n, m, a = map(int, input("Enter the theatre square dimensions (n m) and flagstone size (a), separated by spaces: ").split())

        # Ensure that inputs are positive as per the problem description.
        if not (n > 0 and m > 0 and a > 0):
            print("Error: All dimensions must be positive integers.")
            return

        # 2. Calculate the number of flagstones for the 'n' dimension.
        # This uses the integer arithmetic equivalent of ceiling division,
        # which avoids branching (if/else statements).
        # It's the optimized logic from Question 2 & 3.
        flagstones_n = (n + a - 1) // a

        # 3. Calculate the number of flagstones for the 'm' dimension.
        flagstones_m = (m + a - 1) // a

        # 4. The total number of flagstones is the product of the two dimensions.
        total_flagstones = flagstones_n * flagstones_m
        
        # 5. Output the result showing each number in the final equation.
        print(f"Calculation: {flagstones_n} * {flagstones_m} = {total_flagstones}")
        print(f"Total flagstones needed: {total_flagstones}")

    except ValueError:
        print("Invalid input. Please enter three integers separated by spaces.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Execute the function to solve the problem
solve_theatre_square()
