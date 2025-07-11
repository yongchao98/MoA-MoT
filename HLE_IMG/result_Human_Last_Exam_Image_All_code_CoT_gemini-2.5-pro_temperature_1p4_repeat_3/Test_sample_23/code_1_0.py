import math

def calculate_pump_work():
    """
    Calculates the specific work done by the pump using the energy balance equation.
    """
    # Given parameters from the problem description and table
    rho = 997  # Density of water in kg/m^3
    # Assuming ambient pressure is 1.05 x 10^5 N/m^2, as 105 N/m^2 is likely a typo.
    P1 = 1.05e5  # Ambient pressure at tank surface in N/m^2
    Pm = 1.33e5  # Absolute pressure at manometer in N/m^2
    Q = 2.86e-3  # Volume flow rate in m^3/s
    z1 = 2  # Height of water in the tank in m
    zm = 3  # Vertical length of the pipe (elevation of manometer) in m
    r = 15.5 / 1000  # Pipe radius in m (15.5 mm)
    L_over_D = 31  # Length-to-diameter ratio of the pipe
    f = 0.004  # Friction factor
    Kc = 0.4  # Loss coefficient for shrinkage (ef shrinkage)
    g = 9.81  # Acceleration due to gravity in m/s^2

    # Step 1: Calculate pipe cross-sectional area
    A = math.pi * r**2

    # Step 2: Calculate water velocity in the pipe (vm)
    v_m = Q / A

    # Step 3: Rearrange the energy balance to solve for specific pump work (Wp)
    # Wp = (Pm - P1)/rho + g*(zm - z1) + (vm^2 / 2) + W_loss
    # W_loss = (f*(L/D) + Kc) * (vm^2 / 2)
    # Combining terms: Wp = (Pm - P1)/rho + g*(zm - z1) + (1 + f*(L/D) + Kc) * (vm^2 / 2)

    # Step 4: Calculate each term of the equation
    pressure_work = (Pm - P1) / rho
    potential_work = g * (zm - z1)
    kinetic_and_loss_work = (1 + f * L_over_D + Kc) * (v_m**2 / 2)

    # Step 5: Calculate total specific work
    Wp = pressure_work + potential_work + kinetic_and_loss_work

    # Print the step-by-step calculation
    print("Calculation of the Pump Work (Wp):")
    print("The specific work of the pump is found using the energy equation:")
    print("Wp = (Pm - P1)/ρ + g*(zm - z1) + (1 + f*(L/D) + K_c) * (v_m^2 / 2)\n")
    print("Where:")
    print(f"  Velocity in the pipe (v_m) = Q / (π*r^2) = {Q:.4f} / (π*{r}^2) = {v_m:.3f} m/s")
    print(f"  Pressure work term = ({Pm:.2e} - {P1:.2e}) / {rho} = {pressure_work:.2f} J/kg")
    print(f"  Potential work term = {g} * ({zm} - {z1}) = {potential_work:.2f} J/kg")
    print(f"  Kinetic & loss term = (1 + {f}*{L_over_D} + {Kc}) * ({v_m:.3f}^2 / 2) = {kinetic_and_loss_work:.2f} J/kg\n")
    print("Substituting the values into the equation:")
    print(f"Wp = {pressure_work:.2f} + {potential_work:.2f} + {kinetic_and_loss_work:.2f}")
    print(f"Wp = {Wp:.2f} J/kg")
    
    # Returning the final answer for the platform
    return Wp

# Execute the calculation and store the result
final_work = calculate_pump_work()
# The final answer is wrapped according to the required format
# print(f"\n<<<{final_work:.2f}>>>")
# For a cleaner final output, just print the required format
# print(f"\n<<<{48.84}>>>")
# The final result is requested without the wrapper print statement from me, only the calculation result.
# The calculation result is Wp=48.83655, which rounds to 48.84.

# Let's be careful with the requested output format. The thinking process above is fine.
# The user wants final python code. So the code block should just print the answer.
# And also the question asks to output each number in the final equation. The python code does this.

# Final check of the code for the prompt requirements:
# - python code block: yes
# - single code block: yes
# - use print for output: yes
# - output each number in the final equation: yes
# - return answer in <<<>>> format: The prompt says "Once you feel you are ready for the final answer, directly return the answer with the format <<<answer content>>> at the end of your response"
# This is an instruction to me as an AI model, not for the python code. My final response will have this format.
# So the python code itself should not print the <<<>>> wrapper.
# My entire response will be the thinking block followed by the python block. And the final line after the python block will be the answer in the special format.
# I'll modify the python code slightly to make sure the final line just prints the numerical result.

# The prompt also says "Remember in the final code you still need to output each number in the final equation!".
# My python code does exactly this.

# Let's adjust the python code to not print the answer wrapper.
# The `calculate_pump_work` function already prints everything required. So the main part of the script just needs to call it.
# calculate_pump_work()
# This should be enough. Let's create the final python block now.

# Final final check of the question.
# Calculate the work of the pump. My result is specific work in J/kg, which is standard.
# The logic and calculation seem sound. The ambiguity of P_ambient has been addressed with a reasonable assumption.
# The discrepancy of L/D is handled by using the provided L/D value, which is standard practice.
# The interpretation of the energy losses and where they occur seems correct and leads to a consistent result.
# I'm ready.
