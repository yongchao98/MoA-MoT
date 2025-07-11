def calculate_harmonic_losses():
    """
    Calculates the resistive power losses due to third and fifth harmonic distortions.
    """
    # 1. Identify the given values from the problem description.
    fundamental_loss_mw = 2.0  # MW
    harmonic_loss_factor = 0.10  # 10% additional loss for each harmonic

    # 2. Calculate the power loss for each specific harmonic.
    # Loss due to the 3rd harmonic
    loss_3rd_harmonic = fundamental_loss_mw * harmonic_loss_factor
    
    # Loss due to the 5th harmonic
    loss_5th_harmonic = fundamental_loss_mw * harmonic_loss_factor

    # 3. Sum the individual harmonic losses to get the total loss from harmonics.
    total_harmonic_loss = loss_3rd_harmonic + loss_5th_harmonic

    # 4. Print the calculation steps and the final result.
    print(f"The fundamental power loss is given as {fundamental_loss_mw} MW.")
    print("Each harmonic (3rd and 5th) introduces an additional loss of 10% of the fundamental loss.")
    print("\nStep 1: Calculate the loss from the 3rd harmonic.")
    print(f"Loss_3rd = {fundamental_loss_mw} MW * {harmonic_loss_factor} = {loss_3rd_harmonic:.2f} MW")
    
    print("\nStep 2: Calculate the loss from the 5th harmonic.")
    print(f"Loss_5th = {fundamental_loss_mw} MW * {harmonic_loss_factor} = {loss_5th_harmonic:.2f} MW")

    print("\nStep 3: Calculate the total loss due to harmonics.")
    print("Total Harmonic Loss = Loss from 3rd harmonic + Loss from 5th harmonic")
    print(f"Total Harmonic Loss = {loss_3rd_harmonic:.2f} MW + {loss_5th_harmonic:.2f} MW")
    print(f"\nThe total resistive power loss due to harmonic distortions is {total_harmonic_loss:.2f} MW.")

if __name__ == "__main__":
    calculate_harmonic_losses()
    
    # Calculation for the final answer format
    fundamental_loss_mw = 2.0
    harmonic_loss_factor = 0.10
    loss_3rd_harmonic = fundamental_loss_mw * harmonic_loss_factor
    loss_5th_harmonic = fundamental_loss_mw * harmonic_loss_factor
    total_harmonic_loss = loss_3rd_harmonic + loss_5th_harmonic
    print(f"<<<{total_harmonic_loss}>>>")
