import math

def solve_chain_problem():
    """
    This function calculates the time when the chain first loses contact with the ground.
    It solves the derived quadratic equation for sin(t) and then finds t.
    """
    
    # The problem simplifies to solving a quadratic equation for s = sin(t):
    # a*s^2 + b*s + c = 0
    # where the equation is s^2 - (40 - sqrt(3))s + (39 - 4*sqrt(3)) = 0
    
    # Define the coefficients of the quadratic equation
    a = 1.0
    b = -(40 - math.sqrt(3))
    c = 39 - 4 * math.sqrt(3)
    
    print("The problem reduces to solving a quadratic equation for s = sin(t):")
    print("a*s^2 + b*s + c = 0")
    print("The symbolic equation is: s^2 - (40 - sqrt(3))s + (39 - 4*sqrt(3)) = 0")
    print("\nThe numerical values for the coefficients are:")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")
    
    # Calculate the discriminant
    discriminant = b**2 - 4 * a * c
    
    if discriminant < 0:
        print("\nThe equation has no real solutions for sin(t).")
        return
        
    # Calculate the two possible solutions for s = sin(t)
    s1 = (-b + math.sqrt(discriminant)) / (2 * a)
    s2 = (-b - math.sqrt(discriminant)) / (2 * a)
    
    print(f"\nThe two potential solutions for s = sin(t) are: {s1:.6f} and {s2:.6f}")
    
    # The value of sin(t) must be between -1 and 1. We choose the valid solution.
    valid_s = None
    if -1 <= s1 <= 1:
        valid_s = s1
    elif -1 <= s2 <= 1:
        valid_s = s2
        
    if valid_s is None:
        print("Neither solution for sin(t) is physically possible.")
        return
        
    print(f"\nThe physically valid solution is sin(t) = {valid_s:.6f}")
    
    # Calculate the time 't' using the arcsin function.
    # We take the principal value since we are looking for the first time this occurs.
    time = math.asin(valid_s)
    
    print("\nThe time when the chain first loses contact with the ground is:")
    print(f"{time:.4f} seconds")
    
    # Final answer in the required format
    global final_answer
    final_answer = time

# Execute the solver
final_answer = 0
solve_chain_problem()
print(f'<<<{final_answer:.4f}>>>')