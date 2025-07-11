import math

def calculate_blowup_condition(x0):
    """
    Calculates the critical threshold for y(0) for a given x(0) > 1.
    The solution to the system of ODEs blows up if y(0) is less than this threshold.
    The function also prints the step-by-step calculation.
    """
    if x0 <= 1:
        print("This analysis assumes x(0) > 1.")
        return

    # Let's show the equation first
    print(f"The solution blows up for y(0) < sqrt(2*x(0) - 3*x(0)^(2/3) + 1)")
    print("-" * 20)
    
    # Step-by-step calculation
    print(f"For x(0) = {x0}:")
    
    # 1. Calculate x0^(2/3)
    x0_pow_2_3 = x0**(2/3)
    # The term x0**(1/3) is also calculated implicitly, let's show it for clarity
    x0_pow_1_3 = x0**(1/3)
    print(f"First, calculate the terms involving x(0):")
    print(f"x(0)^(1/3) = {x0}^(1/3) = {x0_pow_1_3:.4f}")
    print(f"x(0)^(2/3) = {x0}^(2/3) = {x0_pow_2_3:.4f}")

    # 2. Substitute into the expression under the square root
    term1 = 2 * x0
    term2 = 3 * x0_pow_2_3
    term3 = 1
    print("\nSubstitute these values into the formula: y(0) < sqrt(2*x(0) - 3*x(0)^(2/3) + 1)")
    # Using format to represent the numbers in the equation
    print(f"y(0) < sqrt(2*({x0}) - 3*({x0_pow_2_3:.4f}) + {term3})")
    
    # 3. Perform the arithmetic
    inside_sqrt = term1 - term2 + term3
    print(f"y(0) < sqrt({term1} - {term2:.4f} + {term3})")
    print(f"y(0) < sqrt({inside_sqrt:.4f})")
    
    # 4. Final result
    ycrit = math.sqrt(inside_sqrt)
    print("\nFinal Result:")
    print(f"The solution blows up for y(0) < {ycrit:.4f}")

# We use x(0) = 8 as an example
x_initial = 8.0
calculate_blowup_condition(x_initial)