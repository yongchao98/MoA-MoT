import sympy

def find_t_bounds():
    """
    This function outlines the derivation and calculation for the bounds of t.
    """
    print("Step 1: Define the intervals and the core condition.")
    print("Let the range of a_i be S = [-1, t].")
    print("The range of sums X = a_0 + a_2 is the interval I = [-2, 2t].")
    print("The condition is: for any x in I, 1/x must also be in I.")
    
    print("\nStep 2: Establish the domain for t.")
    print("To avoid division by zero, x cannot be 0, so I must not contain 0.")
    print("This requires 2t < 0, which means t < 0.")
    print("The interval [-1, t] is valid for t >= -1. So, the domain for t is [-1, 0).")
    
    print("\nStep 3: Set up inequalities from the subset condition.")
    print("The image of I = [-2, 2t] under f(x)=1/x is [1/(2t), -1/2].")
    print("The condition [1/(2t), -1/2] subset of [-2, 2t] yields two inequalities.")
    print("  1) -2 <= 1/(2*t)")
    print("  2) -1/2 <= 2*t")
    
    print("\nStep 4: Solve for t.")
    print("Solving ' -1/2 <= 2*t ' gives t >= -1/4.")
    print("Solving ' -2 <= 1/(2*t) ' for t < 0 gives -4*t >= 1, which simplifies to t <= -1/4.")

    print("\nAt the boundary points of these inequalities, we have the final equations for t:")
    # From -2 <= 1/(2*t), boundary is -2 = 1/(2*t) => -4*t = 1
    final_eq_1 = "-4 * t = 1"
    print(f"Final equation from the first inequality boundary: {final_eq_1}")
    print("The numbers in this equation are -4 and 1.")

    # From -1/2 <= 2*t, boundary is -1/2 = 2*t => -1 = 4*t
    final_eq_2 = "-1 = 4 * t"
    print(f"Final equation from the second inequality boundary: {final_eq_2}")
    print("The numbers in this equation are -1 and 4.")

    print("\nStep 5: Determine the bounds.")
    print("The conditions t >= -1/4 and t <= -1/4 must both be satisfied.")
    print("This is only possible if t is exactly -1/4.")
    
    lower_bound = -0.25
    upper_bound = -0.25
    
    print(f"\nThe lower bound for t is {lower_bound}.")
    print(f"The upper bound for t is {upper_bound}.")
    
    print("\nFinal Answer:")
    print(lower_bound, upper_bound)

find_t_bounds()

# The final answer in the required format
# <<< -0.25 -0.25 >>>