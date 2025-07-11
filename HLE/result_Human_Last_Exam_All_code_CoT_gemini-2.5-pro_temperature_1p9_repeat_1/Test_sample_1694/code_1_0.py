import cmath
from fractions import Fraction

def solve_and_sum_first_coordinate(A_diag, B_diag, C_diag, eq_num):
    """
    Solves one of the matrix equations for the first coordinate and sums the solutions.
    """
    print(f"--- Analyzing Equation {eq_num} ---")
    A_11 = A_diag[0]
    # B is a scalar matrix 6*I, so its diagonal elements are 6
    B_11 = 6
    C_11 = C_diag[0]

    # The system of equations is given by: A * X^2 + X^2 * B = C
    # Let Y = X^2. The equation is AY + YB = C.
    # Since all matrices A, B, C are diagonal, Y must also be diagonal.
    # The equation for the (1,1) element is: A_11 * Y_11 + Y_11 * B_11 = C_11
    # (A_11 + B_11) * Y_11 = C_11
    # Y_11 = C_11 / (A_11 + B_11)
    
    Y_11 = Fraction(C_11, A_11 + B_11)
    
    print(f"The equation for the (1,1) element of Y=X_{eq_num}^2 is ({A_11} + {B_11}) * Y_{{{eq_num},11}} = {C_11}")
    print(f"This gives Y_{{{eq_num},11}} = {C_11} / {A_11 + B_11} = {Y_11}")
    
    # Since X is also diagonal, X_11^2 = Y_11.
    # The solutions for X_11 are sqrt(Y_11) and -sqrt(Y_11).
    sol1 = cmath.sqrt(Y_11)
    sol2 = -sol1
    
    print(f"The solutions for the first coordinate of X_{eq_num} are the square roots of {Y_11}:")
    print(f"Solution 1: {sol1}")
    print(f"Solution 2: {sol2}")
    
    # The sum of these solutions is always zero.
    coord_sum = sol1 + sol2
    # The 'final equation' for this part is the summation
    print(f"The sum of these solutions is: {sol1} + ({sol2}) = {coord_sum}")
    return coord_sum

# Data from the first equation
A1_diag = [5, -5]
C1_diag = [Fraction(-53, 12), 0]
sum1 = solve_and_sum_first_coordinate(A1_diag, [6, 6], C1_diag, 1)

# Data from the second equation
A2_diag = [4, -5]
C2_diag = [Fraction(-3, 11), 0]
sum2 = solve_and_sum_first_coordinate(A2_diag, [6, 6], C2_diag, 2)

# Final total sum
print("\n--- Total Sum ---")
total_sum = sum1 + sum2
# The 'final equation' for the whole problem
print("The final equation for the total sum is the sum of the sums from each equation.")
print(f"Total Sum = {sum1} + {sum2} = {total_sum}")

# The numerical answer is the real part of the total sum.
final_answer = total_sum.real
print(f"\nThe sum of the first coordinate of solutions is {final_answer}.")
<<<0>>>