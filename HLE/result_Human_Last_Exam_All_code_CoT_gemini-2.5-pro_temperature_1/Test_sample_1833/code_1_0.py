import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # Given parameters
    m = 1.0  # mass of the ring in kg
    M = 1.0  # mass of the object in kg
    g = 9.8  # acceleration due to gravity in m/s^2
    theta_deg = 60.0  # angle in degrees

    # Convert angle to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)
    s = math.sin(theta_rad)
    c = math.cos(theta_rad)

    print(f"Solving for m = {m} kg, M = {M} kg, theta = {theta_deg} degrees, and g = {g} m/s^2.\n")

    # Step 1: Calculate the term L*omega^2 (L*ω^2), which represents the contribution of
    # angular velocity to acceleration. L is the string length, ω is the angular velocity.
    # This term is derived from the conservation of energy.
    # Formula: L*ω^2 = (2 * g * sin(θ) * (m + M)) / (m + M * cos(θ)^2)
    
    L_omega_sq_numerator = 2 * g * s * (m + M)
    L_omega_sq_denominator = m + M * c**2
    L_omega_sq = L_omega_sq_numerator / L_omega_sq_denominator

    print("Step 1: Calculate L*ω^2")
    print(f"L*ω^2 = (2 * {g} * sin({theta_deg}) * ({m} + {M})) / ({m} + {M} * cos({theta_deg})^2)")
    print(f"L*ω^2 = ({L_omega_sq_numerator:.4f}) / ({L_omega_sq_denominator:.4f}) = {L_omega_sq:.4f} m/s^2\n")

    # Step 2: Calculate the tension T using the derived formula.
    # This formula relates tension to the forces and accelerations in the system.
    # Formula: T = (m * M / (m + M * cos(θ)^2)) * (g * sin(θ) + L*ω^2)

    term1_numerator = m * M
    term1_denominator = m + M * c**2
    term1 = term1_numerator / term1_denominator
    
    term2 = g * s + L_omega_sq
    
    tension = term1 * term2

    print("Step 2: Calculate the Tension (T)")
    print(f"T = ({m} * {M} / ({m} + {M} * cos({theta_deg})^2)) * ({g} * sin({theta_deg}) + L*ω^2)")
    print(f"T = ({term1_numerator:.2f} / {term1_denominator:.2f}) * ({g * s:.4f} + {L_omega_sq:.4f})")
    print(f"T = {term1:.4f} * {term2:.4f}")
    
    print(f"\nThe final tension is: {tension:.2f} N")
    
    # Returning the final numerical answer as a string for the required format
    return f"{tension:.2f}"

# Execute the function and capture the final answer
final_answer = calculate_tension()
print(f"\n<<<__{final_answer}__>>>")
