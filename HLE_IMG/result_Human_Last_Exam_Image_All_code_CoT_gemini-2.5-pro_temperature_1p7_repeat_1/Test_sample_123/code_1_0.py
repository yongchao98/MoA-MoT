import sympy

def solve_bran_castle_escape():
    """
    Sets up and solves the system of linear equations for the Bran Castle escape problem.
    """
    # Define the 8 unknown probabilities as symbolic variables
    p_me, p_gh, p_ck, p_sp, p_l, p_kh, p_qr, p_tc = sympy.symbols('p_me p_gh p_ck p_sp p_l p_kh p_qr p_tc')
    
    # Define the system of 8 linear equations.
    # Each equation is written in the form f(p_i) = 0.
    # Note: P_TR = 1, P_VL = 0
    
    # 1. Main Entrance
    eq1 = p_me - sympy.Rational(1, 2) * p_gh - sympy.Rational(1, 2) * p_ck
    
    # 2. Great Hall
    eq2 = p_gh - sympy.Rational(1, 4) * p_me - sympy.Rational(1, 4) * p_sp - sympy.Rational(1, 4) * p_l - sympy.Rational(1, 4) * p_kh
    
    # 3. Castle Kitchen
    eq3 = p_ck - sympy.Rational(1, 2) * p_me - sympy.Rational(1, 2) * p_kh
    
    # 4. Secret Passage
    eq4 = p_sp - sympy.Rational(1, 4) * p_gh - sympy.Rational(1, 4) * p_tc - sympy.Rational(1, 4) * 1
    
    # 5. Library
    eq5 = p_l - sympy.Rational(1, 3) * p_gh - sympy.Rational(1, 3) * p_kh - sympy.Rational(1, 3) * p_qr
    
    # 6. Knights' Hall
    eq6 = p_kh - sympy.Rational(1, 4) * p_gh - sympy.Rational(1, 4) * p_l - sympy.Rational(1, 4) * p_qr - sympy.Rational(1, 4) * p_ck
    
    # 7. Queen's Room
    eq7 = p_qr - sympy.Rational(1, 3) * p_l - sympy.Rational(1, 3) * p_kh - sympy.Rational(1, 3) * 1
    
    # 8. Torture Chamber
    eq8 = p_tc - sympy.Rational(1, 2) * p_sp

    # Solve the system of equations
    variables = [p_me, p_gh, p_ck, p_sp, p_l, p_kh, p_qr, p_tc]
    solution = sympy.solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8], variables, dict=True)[0]
    
    p_me_sol = solution[p_me]
    p_gh_sol = solution[p_gh]
    p_ck_sol = solution[p_ck]

    # Print the equation for the Main Entrance and substitute the solved values
    print("The governing equation for the Main Entrance is:")
    print(f"P_ME = (1/2) * P_GH + (1/2) * P_CK\n")
    print("Substituting the solved probabilities:")
    # The requirement is to output each number in the final equation.
    print(f"P_ME = (1 / 2) * ({p_gh_sol}) + (1 / 2) * ({p_ck_sol})")
    
    # Calculate the right hand side to show the equation holds
    rhs = sympy.Rational(1, 2) * p_gh_sol + sympy.Rational(1, 2) * p_ck_sol
    print(f"{p_me_sol} = {rhs}\n")
    
    # Final answer
    print(f"The exact probability of reaching the Treasure Room from the Main Entrance is {p_me_sol}.")

if __name__ == '__main__':
    solve_bran_castle_escape()