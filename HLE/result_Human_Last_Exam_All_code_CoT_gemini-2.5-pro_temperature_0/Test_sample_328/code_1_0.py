import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
output_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer

print("This script determines the best average age estimate for two mating insect species based on their life cycles.")
print("\n--- Analysis of Tridactylophagus tartari (Strepsiptera) ---")
print("1. The adult male of this species is a twisted-wing parasite with an extremely short lifespan.")
print("2. After eclosion (emerging as an adult), its only function is to find a mate and reproduce.")
print("3. It does not feed and dies within hours, almost always in less than a day.")
print("4. Therefore, a reasonable age estimate is a very short period, such as a few hours to a day.")

print("\n--- Analysis of Periplaneta americana (American Cockroach) ---")
print("1. The American cockroach has a much longer lifespan.")
print("2. After eclosion (the final molt to adulthood), it requires a maturation period of about one week to become sexually active.")
print("3. Mating can only occur after this period, so its age must be greater than one week.")
print("4. Considering the time to mature and find a mate, 'one month' is a reasonable average estimate from the choices provided.")

print("\n--- Conclusion ---")
print("We need a pair of estimates: (a few hours to 1 day, ~1 month).")
print("Comparing the options, the pair '16 hours' for Tridactylophagus tartari and 'one month' for Periplaneta americana is the most fitting.")
print("\nFinal Equation:")
# The user requested the final equation with each number
age_tridactylophagus = "16 hours"
age_periplaneta = "one month"
print(f"Estimated age of Tridactylophagus tartari male + Estimated age of Periplaneta americana male = {age_tridactylophagus} + {age_periplaneta}")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
captured_output = output_buffer.getvalue()

# Print the captured output to the actual console
print(captured_output)
print("<<<I>>>")