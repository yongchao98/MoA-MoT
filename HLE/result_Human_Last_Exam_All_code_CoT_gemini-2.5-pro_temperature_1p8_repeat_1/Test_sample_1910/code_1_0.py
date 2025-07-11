import math

def calculate_q_space_peak_location():
    """
    Calculates the Q-space location for the second major diffraction peak of
    orthorhombic NaMgH3.
    """
    # Step 1: Define lattice parameters for orthorhombic NaMgH3 in Angstroms.
    a = 5.755  # Angstrom
    b = 8.046  # Angstrom
    c = 5.688  # Angstrom

    # Step 2: Define Miller indices (h, k, l) for the second major diffraction peak.
    # For NaMgH3 (Pnma), the second most intense peak is the (200) reflection.
    h, k, l = 2, 0, 0

    print("The crystal structure of NaMgH3 is orthorhombic (Pnma).")
    print(f"Lattice parameters: a = {a} Å, b = {b} Å, c = {c} Å.")
    print(f"The second major diffraction peak corresponds to the ({h}{k}{l}) plane.\n")
    
    # Step 3: Calculate the Q-space location.
    # The formula is Q = 2 * pi / d, where 1/d^2 = (h/a)^2 + (k/b)^2 + (l/c)^2.
    # This simplifies to Q = 2 * pi * sqrt((h/a)^2 + (k/b)^2 + (l/c)^2).

    # Calculate the term inside the square root
    d_sq_inv = (h / a)**2 + (k / b)**2 + (l / c)**2
    
    # Calculate Q
    q_value = 2 * math.pi * math.sqrt(d_sq_inv)
    
    # Step 4: Print the final equation and the result.
    print("The Q-space location is calculated using the formula:")
    print("Q = 2 * pi * sqrt((h/a)^2 + (k/b)^2 + (l/c)^2)\n")
    print("Substituting the values:")
    # We explicitly show the numbers in the final equation string.
    equation = f"Q = 2 * pi * sqrt(({h}/{a})^2 + ({k}/{b})^2 + ({l}/{c})^2)"
    result = f"= {q_value:.4f} 1/Angstrom"
    print(equation + " " + result)

# Run the calculation
calculate_q_space_peak_location()

# The final answer for the bot to extract.
# This part will not be printed to the user.
final_answer_value = 2 * math.pi * math.sqrt((2/5.755)**2)
# final_answer_value is approx 2.1836
# <<<2.184>>> # rounded to 3 decimal places
<<<2.1836>>>