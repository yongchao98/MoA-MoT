import cmath

def solve_and_sum():
    """
    Solves the two matrix equations and finds the sum of the first coordinates of the solutions.
    """
    # --- Equation 1 ---
    # The equation is of the form A1*X1^2 + X1^2*B1 = C1
    # Since B1 is 6*I, this simplifies to (A1 + 6*I)*X1^2 = C1
    # Let Y1 = X1^2
    
    # (A1 + 6I)[0,0] is 5 + 6 = 11
    # C1[0,0] is -53/12
    y1_11 = (-53/12) / 11
    
    print("Step 1: Analyzing the first equation.")
    print("The first equation simplifies to (A1 + 6I) * X1^2 = C1.")
    print(f"From this, we find the (1,1) element of X1^2 to be -53/12 / 11 = {y1_11:.4f}")
    
    # The first coordinate of X1, let's call it x1, satisfies x1^2 = y1_11
    sol1_pos = cmath.sqrt(y1_11)
    sol1_neg = -sol1_pos
    sum1 = sol1_pos + sol1_neg
    
    print(f"The possible values for the first coordinate of X1 are the square roots of {y1_11:.4f}, which are:")
    print(f"Solution 1a: {sol1_pos}")
    print(f"Solution 1b: {sol1_neg}")
    print(f"The sum of these two values is {sum1}.\n")

    # --- Equation 2 ---
    # The equation is (A2 + 6*I)*X2^2 = C2
    # Let Y2 = X2^2
    
    # (A2 + 6I)[0,0] is 4 + 6 = 10
    # C2[0,0] is -3/11
    y2_11 = (-3/11) / 10
    
    print("Step 2: Analyzing the second equation.")
    print("The second equation simplifies to (A2 + 6I) * X2^2 = C2.")
    print(f"From this, we find the (1,1) element of X2^2 to be -3/11 / 10 = {y2_11:.4f}")
    
    # The first coordinate of X2, let's call it x2, satisfies x2^2 = y2_11
    sol2_pos = cmath.sqrt(y2_11)
    sol2_neg = -sol2_pos
    sum2 = sol2_pos + sol2_neg

    print(f"The possible values for the first coordinate of X2 are the square roots of {y2_11:.4f}, which are:")
    print(f"Solution 2a: {sol2_pos}")
    print(f"Solution 2b: {sol2_neg}")
    print(f"The sum of these two values is {sum2}.\n")

    # --- Final Sum ---
    total_sum = sum1 + sum2
    print("Step 3: Calculating the total sum.")
    print("The total sum is the sum of all possible first coordinates from both equations.")
    print("The final equation representing the sum is:")
    print(f"Sum = ({sol1_pos}) + ({sol1_neg}) + ({sol2_pos}) + ({sol2_neg})")
    # Clean up floating point representation for the final output
    clean_total_sum = round(total_sum.real) if total_sum.imag == 0 else total_sum
    print(f"Sum = {sum1} + {sum2} = {clean_total_sum}")

    return clean_total_sum

if __name__ == '__main__':
    final_answer = solve_and_sum()
    print(f"\nThe sum of the first coordinate of solutions is {final_answer}.")
    print(f'<<<{final_answer}>>>')