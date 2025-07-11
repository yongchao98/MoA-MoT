import math

def solve_for_theta():
    """
    This function solves for the value of lim <b_p, z_p>.
    """
    # Given values
    # a = lim <h_p, b_p>
    a = 0.9375
    # b = lim <h_p, z_p>
    b = 0.9

    # We derived the quadratic equation for theta = lim <b_p, z_p>:
    # C2*theta^2 + C1*theta + C0 = 0
    # where C2=1, C1 = -2ab, C0 = a^2 + b^2 - 1
    C2 = 1
    C1 = -2 * a * b
    C0 = a**2 + b**2 - 1

    print(f"The values used are:")
    print(f"a = <h, b> = {a}")
    print(f"b = <h, z> = {b}")
    print("\nThis leads to the quadratic equation for theta = <b, z>:")
    print(f"{C2} * theta^2 + ({C1}) * theta + ({C0}) = 0")

    # Calculate the discriminant
    discriminant = C1**2 - 4 * C2 * C0

    if discriminant < 0:
        print("\nThe equation has no real solutions.")
        return

    # Solve for the two roots
    sol1 = (-C1 + math.sqrt(discriminant)) / (2 * C2)
    sol2 = (-C1 - math.sqrt(discriminant)) / (2 * C2)

    print(f"\nThe two possible solutions for theta are:")
    print(f"Solution 1: {sol1}")
    print(f"Solution 2: {sol2}")
    
    # To choose the correct solution, we enforce the condition that the coefficients
    # alpha and gamma are positive. This requires theta < b/a.
    selection_threshold = b / a
    print(f"\nTo select the correct root, we use the condition theta < b/a = {selection_threshold:.4f}")

    final_solution = None
    if sol1 < selection_threshold:
        final_solution = sol1
        print(f"Solution 1 ({sol1:.4f}) satisfies the condition.")
    
    if sol2 < selection_threshold:
        # Check if we already found a solution. This would be unexpected.
        if final_solution is not None:
             print("Warning: Both solutions satisfy the condition. This is unexpected.")
        final_solution = sol2
        print(f"Solution 2 ({sol2:.4f}) satisfies the condition.")

    if final_solution is not None:
        print(f"\nThe final answer is {final_solution}")
    else:
        print("\nNeither solution satisfies the condition.")

    return final_solution

# Run the solver
final_answer = solve_for_theta()
