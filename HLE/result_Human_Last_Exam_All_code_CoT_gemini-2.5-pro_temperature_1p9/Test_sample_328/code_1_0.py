import textwrap

def solve_mating_age_puzzle():
    """
    Analyzes the mating ages of two insect species to find the best-fit answer
    from a list of choices.
    """

    # Scientific facts about the two species' mating habits post-eclosion (emergence as adult).
    tridactylophagus_info = {
        "name": "Tridactylophagus tartari (a twisted-wing parasite)",
        "fact": "Adult males have an extremely short lifespan, typically lasting only a few hours (e.g., less than 6). Their sole purpose is to mate immediately after emerging."
    }

    periplaneta_info = {
        "name": "Periplaneta americana (American cockroach)",
        "fact": "After molting into an adult, males become sexually mature in a few days. Mating is possible as early as two days post-eclosion, with frequency increasing thereafter."
    }

    # The chosen answer based on the biological facts.
    best_choice_label = "K"
    best_choice_values = ("6 hours", "two days")

    # --- Print the step-by-step reasoning ---
    print("Step 1: Analyzing the Mating Age of Tridactylophagus tartari")
    print("-" * 60)
    print(textwrap.fill(f"The first species is {tridactylophagus_info['name']}. {tridactylophagus_info['fact']}", 70))
    print("Conclusion for species 1: The best estimate for its age must be on the order of hours.\n")


    print("Step 2: Analyzing the Mating Age of Periplaneta americana")
    print("-" * 60)
    print(textwrap.fill(f"The second species is {periplaneta_info['name']}. {periplaneta_info['fact']}", 70))
    print("Conclusion for species 2: The best estimate for its age is on the order of a few days.\n")

    print("Step 3: Selecting the Best Answer")
    print("-" * 60)
    print("We must find the answer choice that has a value in 'hours' for the first species and 'days' for the second.")
    print("Choice K is ('6 hours', 'two days').")
    print("  - '6 hours' is a perfect estimate for Tridactylophagus tartari.")
    print("  - 'two days' is a perfect estimate for Periplaneta americana mating shortly after eclosion.\n")

    print("Final Answer Derivation:")
    print("-" * 60)
    print(f"The best average estimate for the male {tridactylophagus_info['name'].split(' ')[0]}'s age is {best_choice_values[0]}.")
    print(f"The best average estimate for the male {periplaneta_info['name'].split(' ')[0]}'s age is {best_choice_values[1]}.")

if __name__ == "__main__":
    solve_mating_age_puzzle()
    print("\n<<<K>>>")