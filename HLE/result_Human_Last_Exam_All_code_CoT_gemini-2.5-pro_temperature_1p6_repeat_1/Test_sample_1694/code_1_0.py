import cmath

def solve_system():
    """
    Solves the given system of matrix equations and finds the sum of the first coordinates of the solutions.
    """
    
    # --- Equation 1 ---
    # (A1 + 6I) * X1^2 = C1
    # Let Y1 = X1^2. We have [[11, 0], [0, 1]] * Y1 = [[-53/12, 0], [0, 0]]
    # This implies Y1 = [[(-53/12)/11, 0], [0, 0]]
    
    y1_11 = -53.0 / (12.0 * 11.0)
    
    # From X1^2 = Y1, and since Y1 is diagonal, we can find the first coordinate of X1.
    # Let the first coordinate of X1 be x1_11. Then x1_11^2 = y1_11.
    # The two possible solutions for x1_11 are:
    x1_coord1 = cmath.sqrt(y1_11)
    x1_coord2 = -x1_coord1

    # --- Equation 2 ---
    # (A2 + 6I) * X2^2 = C2
    # Let Y2 = X2^2. We have [[10, 0], [0, 1]] * Y2 = [[-3/11, 0], [0, 0]]
    # This implies Y2 = [[(-3/11)/10, 0], [0, 0]]

    y2_11 = -3.0 / (11.0 * 10.0)

    # From X2^2 = Y2, we find the first coordinate of X2.
    # Let the first coordinate of X2 be x2_11. Then x2_11^2 = y2_11.
    # The two possible solutions for x2_11 are:
    x2_coord1 = cmath.sqrt(y2_11)
    x2_coord2 = -x2_coord1

    # The four possible first coordinates from the solutions are:
    solutions = [x1_coord1, x1_coord2, x2_coord1, x2_coord2]

    # Calculate the sum
    total_sum = sum(solutions)
    
    # Print the equation representing the sum of the coordinates
    print("The sum of the first coordinate of the solutions is calculated as follows:")
    # The f-string formatting {c.real:+.4f}{c.imag:+.4f}j shows complex numbers nicely.
    print(f"({x1_coord1.real:+.4f}{x1_coord1.imag:+.4f}j) + ({x1_coord2.real:+.4f}{x1_coord2.imag:+.4f}j) + ({x2_coord1.real:+.4f}{x2_coord1.imag:+.4f}j) + ({x2_coord2.real:+.4f}{x2_coord2.imag:+.4f}j) = {total_sum.real:+.4f}{total_sum.imag:+.4f}j")
    
    return total_sum

result = solve_system()

# Final answer formatting
# Check if the imaginary part is close to zero before casting to float
final_answer = result.real if abs(result.imag) < 1e-9 else result
print(f"\nThe final answer is {final_answer}")
print(f"<<<{final_answer}>>>")