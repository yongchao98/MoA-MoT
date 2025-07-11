import fractions

# --- Step 1: Find the total number of pipes ---

# The number of pipes that are still in tune.
in_tune_count = 200

# The fractions of pipes that fell out of tune.
fraction_fell = fractions.Fraction(1, 3)
fraction_moon = fractions.Fraction(2, 5)

# Calculate the total fraction of pipes that are out of tune.
total_out_of_tune_fraction = fraction_fell + fraction_moon

# The remaining fraction represents the pipes that are still in tune.
in_tune_fraction = 1 - total_out_of_tune_fraction

# Calculate the total number of pipes using the count of in-tune pipes.
# (in_tune_fraction * total_pipes) = in_tune_count
# So, total_pipes = in_tune_count / in_tune_fraction
total_pipes = int(in_tune_count / in_tune_fraction)


# --- Step 2: Calculate the "Lost" Pipes ---

# The number of "lost" (out-of-tune) pipes is the total minus the in-tune pipes.
out_of_tune_count = total_pipes - in_tune_count


# --- Step 3 & 4: Find the final answer and display the equation ---

# The tuner must find half of the lost pipes.
pipes_to_tune = int(out_of_tune_count / 2)

print(f"First, we found the total number of pipes to be {total_pipes}.")
print(f"From that, the number of out-of-tune pipes is {total_pipes} - {in_tune_count} = {out_of_tune_count}.")
print("\nThe question asks how many pipes the tuner must find, which is half of the lost pipes.")
print("The final equation is:")
print(f"{out_of_tune_count} / 2 = {pipes_to_tune}")

print(f"\nThus, the tuner must find {pipes_to_tune} pipes.")
<<<275>>>