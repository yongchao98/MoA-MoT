#
# A python script to calculate the total resistive power loss including harmonics.
#

# 1. Define the given values based on the problem description.
# The fundamental power loss is given as 2 MW.
fundamental_loss = 2.0  # in MW

# The additional power loss for each harmonic order is 10%.
harmonic_loss_percentage = 0.10

# 2. Calculate the power loss for each specified harmonic.
# The problem asks for the third and fifth harmonics.
loss_3rd_harmonic = fundamental_loss * harmonic_loss_percentage
loss_5th_harmonic = fundamental_loss * harmonic_loss_percentage

# 3. Calculate the total power loss.
# This is the sum of the fundamental loss and the losses from the harmonics.
total_loss = fundamental_loss + loss_3rd_harmonic + loss_5th_harmonic

# 4. Print the final equation and the result.
# The output clearly shows each component of the total loss.
print("Calculating the total resistive power loss:")
print(f"Total Loss = Fundamental Loss + 3rd Harmonic Loss + 5th Harmonic Loss")
print(f"Total Loss = {fundamental_loss:.1f} MW + {loss_3rd_harmonic:.1f} MW + {loss_5th_harmonic:.1f} MW")
print(f"Total Loss = {total_loss:.1f} MW")