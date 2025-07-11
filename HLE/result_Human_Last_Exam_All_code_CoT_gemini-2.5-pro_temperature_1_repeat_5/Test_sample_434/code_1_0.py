import math

def solve_blowup_condition():
    """
    Calculates the critical value for y(0) that determines whether the solution
    to the given system of differential equations will blow up, for a specific
    value of x(0) > 1.
    """
    # Assume a value for x(0) > 1. We choose x(0) = 8 for a clean cube root.
    x0 = 8.0

    print(f"Analyzing the system for the initial condition x(0) = {x0}.")
    print("The condition for blow-up depends on a critical threshold for y(0).")
    print("This threshold is given by the separatrix curve passing through the saddle point (1,0).")
    print("The equation for the stable branch of the separatrix for x > 1 is:")
    print("y_crit = (x^(1/3) - 1) * sqrt(2*x^(1/3) + 1)\n")

    # Calculate the components of the equation
    z = x0**(1/3)
    term1 = z - 1
    term2_inner = 2 * z + 1
    term2 = math.sqrt(term2_inner)
    y_crit = term1 * term2

    # Print the step-by-step calculation
    print(f"For x(0) = {x0}:")
    print(f"  x(0)^(1/3) = {x0}^(1/3) = {z}")
    print(f"  y_crit = ({z} - 1) * sqrt(2 * {z} + 1)")
    print(f"         = {term1} * sqrt({term2_inner})")
    print(f"         = {term1} * {term2}")
    print(f"         = {y_crit}\n")

    # Print the final condition
    print("The solution blows up if y(0) is less than this critical value.")
    print(f"Blow-up condition for x(0) = {x0}: y(0) < {y_crit}")

solve_blowup_condition()

# The final answer is the numerical threshold for y(0) when x(0) = 8.
# y_crit = (8^(1/3)-1)*sqrt(2*8^(1/3)+1) = (2-1)*sqrt(2*2+1) = sqrt(5)
final_answer = math.sqrt(5)
print(f"\n<<<y(0) < {final_answer}>>>")