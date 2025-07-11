# Define the given values for the calculation
fundamental_loss_mw = 2.0  # in Megawatts
additional_loss_rate = 0.10  # 10% additional loss for each harmonic

# Calculate the additional loss for each harmonic
# Each harmonic adds 10% of the fundamental loss.
loss_3rd_harmonic = fundamental_loss_mw * additional_loss_rate
loss_5th_harmonic = fundamental_loss_mw * additional_loss_rate

# Calculate the total power loss by summing the fundamental loss and harmonic losses
total_loss_mw = fundamental_loss_mw + loss_3rd_harmonic + loss_5th_harmonic

# Print the step-by-step calculation
print("Calculation of Total Resistive Power Loss:")
print("-" * 45)
print(f"Fundamental Power Loss: {fundamental_loss_mw:.2f} MW")
print(f"Additional Loss from 3rd Harmonic: {loss_3rd_harmonic:.2f} MW")
print(f"Additional Loss from 5th Harmonic: {loss_5th_harmonic:.2f} MW")
print("-" * 45)

# Print the final equation with all its numerical components
print("Total Power Loss = Fundamental Loss + 3rd Harmonic Loss + 5th Harmonic Loss")
print(f"Total Power Loss = {fundamental_loss_mw:.2f} MW + {loss_3rd_harmonic:.2f} MW + {loss_5th_harmonic:.2f} MW")
print(f"Final Total Power Loss = {total_loss_mw:.2f} MW")

<<<2.4>>>