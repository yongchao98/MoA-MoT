import sys

def calculate_harmonic_power_loss():
    """
    Calculates the resistive power losses due to third and fifth harmonic distortions.
    """
    # Given values
    fundamental_loss = 2.0  # in MW
    harmonic_loss_percentage = 10.0 # in percent

    # Convert percentage to a decimal factor
    harmonic_loss_factor = harmonic_loss_percentage / 100.0

    # Calculate the power loss introduced by each harmonic
    # Each is an additional 10% of the fundamental loss
    loss_3rd_harmonic = fundamental_loss * harmonic_loss_factor
    loss_5th_harmonic = fundamental_loss * harmonic_loss_factor

    # The total loss due to harmonics is the sum of the individual harmonic losses
    total_harmonic_loss = loss_3rd_harmonic + loss_5th_harmonic

    print("Step 1: Calculate the power loss due to the 3rd harmonic.")
    print(f"Loss_3rd = Fundamental Loss * Harmonic Loss Factor")
    print(f"Loss_3rd = {fundamental_loss:.1f} MW * {harmonic_loss_factor:.2f} = {loss_3rd_harmonic:.2f} MW\n")

    print("Step 2: Calculate the power loss due to the 5th harmonic.")
    print(f"Loss_5th = Fundamental Loss * Harmonic Loss Factor")
    print(f"Loss_5th = {fundamental_loss:.1f} MW * {harmonic_loss_factor:.2f} = {loss_5th_harmonic:.2f} MW\n")

    print("Step 3: Calculate the total power loss due to the harmonics.")
    print("Total Harmonic Loss = Loss due to 3rd harmonic + Loss due to 5th harmonic")
    print(f"Total Harmonic Loss = {loss_3rd_harmonic:.2f} MW + {loss_5th_harmonic:.2f} MW")
    
    # Final result
    print(f"\nThe total resistive power loss due to third and fifth harmonic distortions is {total_harmonic_loss:.2f} MW.")
    
    # Capturing the final numerical answer for the platform
    # The final answer is the total harmonic loss in MW.
    final_answer = total_harmonic_loss
    # Use a format that is easily parsable
    sys.stdout.write(f"\n<<<{final_answer:.1f}>>>")

if __name__ == "__main__":
    calculate_harmonic_power_loss()