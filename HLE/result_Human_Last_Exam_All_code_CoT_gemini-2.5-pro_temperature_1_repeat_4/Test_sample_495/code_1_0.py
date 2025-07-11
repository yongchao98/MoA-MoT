def calculate_energy_loss():
    """
    Calculates the energy loss per centimeter for alpha-particles in standard air
    based on Geiger's range-energy relationship.
    """
    # Initial conditions
    E0 = 8.5  # Initial energy in MeV
    R0 = 8.3  # Total range in cm
    x = 4.0   # Distance from the source in cm

    print(f"Given initial energy E₀ = {E0} MeV and total range R₀ = {R0} cm.")
    print(f"We need to find the energy loss per cm at a distance x = {x} cm.\n")
    
    # Step 1: Calculate the constant k from R = k * E^(3/2)
    # k = R0 / E0^(3/2)
    k = R0 / (E0 ** 1.5)
    print(f"Step 1: Calculate the constant k.")
    print(f"k = R₀ / E₀^(3/2) = {R0:.1f} / ({E0:.1f}^1.5) = {k:.4f} cm/MeV^(3/2)\n")

    # Step 2: Calculate the remaining range after traveling x cm
    R_rem = R0 - x
    print(f"Step 2: Calculate the remaining range at x = {x:.1f} cm.")
    print(f"R_rem = R₀ - x = {R0:.1f} - {x:.1f} = {R_rem:.1f} cm\n")

    # Step 3: Calculate the energy E at distance x using the remaining range
    # E = (R_rem / k)^(2/3)
    E_at_x = (R_rem / k) ** (2.0 / 3.0)
    print(f"Step 3: Calculate the energy E at x = {x:.1f} cm.")
    print(f"E = (R_rem / k)^(2/3) = ({R_rem:.1f} / {k:.4f})^(2/3) = {E_at_x:.3f} MeV\n")

    # Step 4: Calculate the energy loss per centimeter, dE/dx
    # dE/dx = 1 / ( (3/2) * k * E^(1/2) )
    dE_dx = 1 / (1.5 * k * (E_at_x ** 0.5))
    print(f"Step 4: Calculate the energy loss per cm (dE/dx).")
    print(f"dE/dx = 1 / (1.5 * k * E^0.5) = 1 / (1.5 * {k:.4f} * {E_at_x:.3f}^0.5)")
    print(f"dE/dx = {dE_dx:.3f} MeV/cm\n")

    print(f"The final calculated energy loss per centimetre is {dE_dx:.3f} MeV/cm.")
    
    return dE_dx

# Run the calculation and store the final answer
final_answer = calculate_energy_loss()
# The final answer format is not printed to the console, it is a marker for the system.
# The value is printed above for the user.
final_answer_formatted_string = f"<<<{final_answer:.3f}>>>"

# The following print is just for capturing the value for the final answer block.
# print(final_answer_formatted_string)