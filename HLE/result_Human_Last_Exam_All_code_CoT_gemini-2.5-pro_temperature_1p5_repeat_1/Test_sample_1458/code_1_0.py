# The number of pipes that still sing pure
in_tune_count = 200

# The fractions of pipes that fell out of tune
fraction1 = 1/3
fraction2 = 2/5

# Calculate the total fraction of out-of-tune pipes
total_fraction_out_of_tune = fraction1 + fraction2

# Calculate the fraction of pipes that are still in tune
fraction_in_tune = 1 - total_fraction_out_of_tune

# Calculate the total number of pipes in the cathedral
total_pipes = in_tune_count / fraction_in_tune

# Calculate the number of "lost" (out-of-tune) pipes
lost_pipes_count = total_pipes - in_tune_count

# The tuner must find half of the lost pipes
pipes_to_find = lost_pipes_count / 2

# The information about three-sevenths and one-fourth is extra detail not needed for the final calculation.

# We will now print the final equation step-by-step with all the numbers.
print("Here is the equation to find the number of pipes the tuner must realign:")
print(f"(({in_tune_count} / (1 - ({int(fraction1*15)}/15 + {int(fraction2*15)}/15))) - {in_tune_count}) / 2 = {int(pipes_to_find)}")