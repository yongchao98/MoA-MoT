import fractions

# Step 1: Define the known values from the poem.
pipes_in_tune = 200
# The fractions of pipes that fell out of tune.
fraction_out_1 = fractions.Fraction(1, 3)
fraction_out_2 = fractions.Fraction(2, 5)

# Step 2: Calculate the total fraction of pipes that are out of tune.
total_fraction_out_of_tune = fraction_out_1 + fraction_out_2

# Step 3: Calculate the fraction of pipes that are still in tune.
fraction_in_tune = 1 - total_fraction_out_of_tune

# Step 4: Calculate the total number of pipes in the cathedral.
# total_pipes = pipes_in_tune / fraction_in_tune
total_pipes = int(pipes_in_tune / fraction_in_tune)

# Step 5: Calculate the number of pipes that are out of tune.
pipes_out_of_tune = total_pipes - pipes_in_tune

# Step 6: Calculate the number of pipes the tuner must find.
# This is half of the out-of-tune pipes.
pipes_to_realign = pipes_out_of_tune / 2

# Print the step-by-step solution
print("Here is the solution to the riddle of the Cathedral's Echo:\n")

print(f"1. First, we find the total fraction of pipes that are out of tune:")
print(f"   1/3 + 2/5 = {fraction_out_1.numerator}/{fraction_out_1.denominator} + {fraction_out_2.numerator}/{fraction_out_2.denominator} = {total_fraction_out_of_tune.numerator}/{total_fraction_out_of_tune.denominator}")

print(f"\n2. Next, we find the fraction of pipes that remain in tune:")
print(f"   1 - {total_fraction_out_of_tune.numerator}/{total_fraction_out_of_tune.denominator} = {fraction_in_tune.numerator}/{fraction_in_tune.denominator}")

print(f"\n3. With {pipes_in_tune} pipes still singing pure, we can find the total number of pipes:")
print(f"   Total Pipes = {pipes_in_tune} / ({fraction_in_tune.numerator}/{fraction_in_tune.denominator}) = {total_pipes}")

print(f"\n4. Now, we calculate the number of pipes that are out of tune:")
print(f"   {total_pipes} (Total) - {pipes_in_tune} (In Tune) = {pipes_out_of_tune}")

print(f"\n5. Finally, the tuner must find half of these 'lost' pipes to realign:")
print(f"   {pipes_out_of_tune} / 2 = {int(pipes_to_realign)}")

print("\nThe tuner must find 275 pipes.")
<<<275>>>