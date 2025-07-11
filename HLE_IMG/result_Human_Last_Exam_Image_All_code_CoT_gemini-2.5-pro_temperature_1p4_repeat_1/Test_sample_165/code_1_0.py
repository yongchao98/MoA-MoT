import sys

def calculate_harmonic_losses():
    """
    Calculates the total resistive power losses on a transmission line,
    including the effects of third and fifth harmonics.
    """
    # 1. Define the given values
    fundamental_loss_mw = 2.0
    harmonic_loss_percentage = 0.10  # 10%

    # 2. Calculate the additional loss for the 3rd harmonic
    third_harmonic_loss_mw = fundamental_loss_mw * harmonic_loss_percentage

    # 3. Calculate the additional loss for the 5th harmonic
    fifth_harmonic_loss_mw = fundamental_loss_mw * harmonic_loss_percentage

    # 4. Calculate the total power loss
    total_loss_mw = fundamental_loss_mw + third_harmonic_loss_mw + fifth_harmonic_loss_mw

    # 5. Print the breakdown and the final equation
    print("This script calculates the total power loss including harmonic distortions.")
    print(f"Fundamental Power Loss (P_1) = {fundamental_loss_mw} MW")
    print(f"Additional loss from 3rd harmonic (P_3) = {harmonic_loss_percentage*100:.0f}% of P_1 = {third_harmonic_loss_mw} MW")
    print(f"Additional loss from 5th harmonic (P_5) = {harmonic_loss_percentage*100:.0f}% of P_1 = {fifth_harmonic_loss_mw} MW")
    print("\nThe total power loss is the sum of the fundamental and harmonic losses.")
    print("\nFinal Equation:")
    # Using print to output the equation with all the numbers
    print(f"Total Loss = {fundamental_loss_mw} MW + {third_harmonic_loss_mw} MW + {fifth_harmonic_loss_mw} MW = {total_loss_mw} MW")

    # Output the final numerical answer in the specified format
    sys.stdout.write(f"<<<{total_loss_mw}>>>")

# Execute the function
calculate_harmonic_losses()