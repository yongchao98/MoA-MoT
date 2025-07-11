import math

def find_astroid_arc_parameter():
    """
    Calculates the parameter 'a' for an astroid arc of a given length.
    The arc is defined by x=(cos(t))^3, y=(sin(t))^3, with 0 <= x <= a.
    The given length of the arc is 3/2.
    """
    # Given arc length
    L_given = 3/2

    print("The problem is to find the value of 'a' for the arc of the astroid")
    print("x = (cos(t))^3, y = (sin(t))^3, with 0 <= x <= a,")
    print(f"given that the total arc length is {L_given}.")
    print("-" * 40)

    print("Step 1: State the formula for the arc length L.")
    print("The total length of the arc defined by 0 <= x <= a is given by the formula:")
    print("L = 3 * (cos(t_a))^2, where a = (cos(t_a))^3.")
    print("-" * 40)

    print("Step 2: Set up the equation with the given length.")
    # The equation is 3 * (cos(t_a))^2 = 3/2
    lhs_factor = 3
    print(f"We have the equation: {lhs_factor} * (cos(t_a))^2 = {L_given}")
    print("-" * 40)

    print("Step 3: Solve for cos(t_a).")
    # (cos(t_a))^2 = (3/2) / 3 = 1/2
    cos_ta_sq_val = L_given / lhs_factor
    print(f"(cos(t_a))^2 = {L_given} / {lhs_factor} = {cos_ta_sq_val}")
    
    # cos(t_a) = sqrt(1/2)
    cos_ta_val = math.sqrt(cos_ta_sq_val)
    print(f"Taking the positive square root, cos(t_a) = sqrt({cos_ta_sq_val}) = 1/sqrt(2)")
    print("-" * 40)

    print("Step 4: Solve for 'a'.")
    print("The constant 'a' is related to t_a by a = (cos(t_a))^3.")
    
    # a = (1/sqrt(2))^3 = 1/(2*sqrt(2)) = sqrt(2)/4
    a_val = cos_ta_val**3
    
    print("Substituting the value of cos(t_a) into the equation for 'a':")
    
    # Printing each number in the final equation as requested
    num1 = 1
    den1_base = 2
    num2 = 1
    den2_base = 2
    num3_sqrt = 2
    den3 = 4
    
    print(f"a = ({num1}/sqrt({den1_base}))^3 = {num2}/({den2_base}*sqrt({den2_base})) = sqrt({num3_sqrt})/{den3}")
    print(f"\nThe numerical value of a is approximately {a_val:.5f}")

find_astroid_arc_parameter()