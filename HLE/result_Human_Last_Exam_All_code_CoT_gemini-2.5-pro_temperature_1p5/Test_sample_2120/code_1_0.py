import sympy

def solve_problem():
    """
    Calculates the average value of the complex coordinates z where the inverse
    of the gradient of S(z) approaches zero.
    """
    # Step 1: Identify the poles of E(z)
    # The poles of E(z) are inferred from the singularities of the arguments
    # in its governing equation: E(f(z)) and E(g(z)), where f(z) has a pole
    # at z=2 and g(z) has a pole at z=1.
    poles_E = [1, 2]
    num_poles_E = len(poles_E)
    sum_poles_E = sum(poles_E)

    print(f"The poles from E(z) are located at z = {poles_E[0]} and z = {poles_E[1]}.")
    print(f"The number of poles from E(z) is {num_poles_E}.")
    print(f"The sum of poles from E(z) is {sum_poles_E}.")
    print("-" * 20)

    # Step 2: Identify the poles of B(z)
    # The poles of B(z) are inferred from the roots of the denominator
    # of the right-hand side of its governing equation: 4*z**4 - z**3 + z**2 + 1 = 0.
    # We use Vieta's formulas to find the sum of these roots without solving for them.
    # For a polynomial a*z**4 + b*z**3 + c*z**2 + d*z + e = 0, sum of roots = -b/a.
    a = 4
    b = -1
    num_poles_B = 4
    sum_poles_B = -b / a

    print(f"The poles from B(z) are the roots of the equation {a}z^4 + ({b})z^3 + ... = 0.")
    print(f"The number of poles from B(z) is {num_poles_B}.")
    print(f"Using Vieta's formulas, the sum of poles from B(z) is -({b})/({a}) = {sympy.S(sum_poles_B)}.")
    print("-" * 20)

    # Step 3: Calculate the average of all poles
    total_num_poles = num_poles_E + num_poles_B
    total_sum_poles = sum_poles_E + sum_poles_B
    
    # Represent numbers as fractions for exact arithmetic
    sum_E_frac = sympy.S(sum_poles_E)
    sum_B_frac = sympy.S(sum_poles_B)
    total_sum_frac = sum_E_frac + sum_B_frac

    print(f"The total number of poles for S(z) is {num_poles_E} + {num_poles_B} = {total_num_poles}.")
    print(f"The total sum of all poles is {sum_E_frac} + {sum_B_frac} = {total_sum_frac}.")
    
    average_value = total_sum_frac / total_num_poles
    
    print("-" * 20)
    print("The final calculation for the average value is:")
    print(f"Average = (Sum of all poles) / (Number of all poles)")
    # Final equation with each number explicitly shown
    print(f"Average = ({sum_poles_E} + {sum_poles_B}) / ({num_poles_E} + {num_poles_B}) = {total_sum_poles} / {total_num_poles} = {average_value}")

solve_problem()
print("\n<<<13/24>>>")