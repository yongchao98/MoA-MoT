import math

def solve_bvp():
    """
    This function solves the given boundary-value problem by printing out the step-by-step derivation.
    """
    # Step 1: Define the constants from the problem description
    A = 10**10
    B = 1/100000 - 1
    T = 10**20
    alpha1 = 0
    # The value of alpha2 is not needed for the solution of X_0(t)

    print("Step 1: Isolate and simplify the equation for X_0(t).")
    print("The differential equation for X_0(t) is: X'_0(t) = -(B + 1)X_0(t) + A")
    print(f"The given values are A = {A} and B = 1/100000 - 1 = {B}.")
    B_plus_1 = B + 1
    print(f"First, we calculate the term (B + 1): ({B}) + 1 = {B_plus_1}.")
    print(f"Substituting the values, the ODE becomes: X'_0(t) = -{B_plus_1}*X_0(t) + {A}.")
    print("This is a first-order linear ODE that can be solved independently of Y_0(t).")
    print("-" * 50)

    print("Step 2: Find the general solution for X_0(t).")
    print("The equation can be written as: X'_0(t) + " + str(B_plus_1) + "*X_0(t) = " + str(A))
    print("The general solution is the sum of a particular solution (X_p) and the homogeneous solution (X_h).")
    # Calculate the particular solution
    X_p = A / B_plus_1
    print(f"The particular solution is a constant X_p = A / (B + 1) = {A} / {B_plus_1} = {X_p:.0e}.")
    print(f"The homogeneous solution is X_h(t) = C1 * exp(-{B_plus_1}*t), where C1 is the integration constant.")
    print(f"The general solution is therefore: X_0(t) = {X_p:.0e} + C1 * exp(-{B_plus_1}*t).")
    print("-" * 50)

    print("Step 3: Apply the boundary condition to find the constant C1.")
    print(f"The boundary condition for X_0 is: X_0(0) - X_0(T) = alpha_1, with T = {T:.0e} and alpha_1 = {alpha1}.")
    print(f"Substituting the general solution into the condition:")
    print(f"X_0(0) = {X_p:.0e} + C1 * exp(0) = {X_p:.0e} + C1")
    exponent = -B_plus_1 * T
    print(f"X_0(T) = {X_p:.0e} + C1 * exp(-{B_plus_1} * {T:.0e}) = {X_p:.0e} + C1 * exp({exponent:.0e})")
    print(f"So, ( {X_p:.0e} + C1 ) - ( {X_p:.0e} + C1 * exp({exponent:.0e}) ) = 0.")
    print(f"This simplifies to: C1 * (1 - exp({exponent:.0e})) = 0.")
    print(f"Since T = {T:.0e} is very large, the term exp({exponent:.0e}) is extremely close to 0, so (1 - exp(...)) is not zero.")
    print("Therefore, for the equation to hold, the constant C1 must be 0.")
    print("-" * 50)

    print("Step 4: Determine the final solution and the value at t = 10^20.")
    print("With C1 = 0, the specific solution for X_0(t) is a constant:")
    final_X0_value = int(X_p)
    print(f"X_0(t) = {final_X0_value:.0e}")
    print(f"The question asks for the value of X_0(t) at t = {T:.0e}.")
    print(f"Since X_0(t) is a constant, its value at t = {T:.0e} is the constant value itself.")
    
    # Final output showing the calculation with the original numbers
    A_int = 10**10
    B_plus_1_num = 1
    B_plus_1_den = 100000
    print("\nThe final calculation is:")
    print(f"X_0({int(T)}) = A / (B + 1) = {A_int} / ({B_plus_1_num}/{B_plus_1_den}) = {final_X0_value}")
    
solve_bvp()