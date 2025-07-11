import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

def solve():
    """
    This function encapsulates the logic to determine the estimated age of mating males
    for Tridactylophagus tartari and Periplaneta americana.
    """
    # Step 1: Analyze Tridactylophagus tartari (a parasitoid wasp).
    # Research indicates these wasps have a very short adult lifespan, typically only 3-4 days.
    # Critically, males often emerge 1-2 days before the females from the same host and then wait.
    # Mating occurs shortly after the female emerges. This means the male's age at mating
    # is at least 1-2 days. Given its entire lifespan is only 3-4 days, an average age of
    # three days is a plausible estimate among the choices.
    age_wasp_str = "three days"
    age_wasp_val = 3

    # Step 2: Analyze Periplaneta americana (the American cockroach).
    # This insect has a long adult lifespan, often exceeding a year. After its final molt
    # (eclosion), it requires about a week to become sexually mature. Mating can then occur
    # throughout its long adult life. Therefore, an average age for a randomly found
    # mating male would likely be several months. Six months is a reasonable estimate.
    age_cockroach_str = "six months"
    age_cockroach_val = 6 # in months

    # Step 3: Print the findings.
    print(f"Based on their life cycles, the best average estimate of the male's age since eclosion is:")
    print(f"- For Tridactylophagus tartari: {age_wasp_str}")
    print(f"- For Periplaneta americana: {age_cockroach_str}")

    # Step 4: Display the final numbers as requested.
    print("\nFinal Answer Breakdown:")
    print(f"Tridactylophagus tartari age estimate (days) = {age_wasp_val}")
    print(f"Periplaneta americana age estimate (months) = {age_cockroach_val}")

solve()

# Restore original stdout
sys.stdout = original_stdout
# Get the output from the buffer
output = captured_output.getvalue()
print(output)