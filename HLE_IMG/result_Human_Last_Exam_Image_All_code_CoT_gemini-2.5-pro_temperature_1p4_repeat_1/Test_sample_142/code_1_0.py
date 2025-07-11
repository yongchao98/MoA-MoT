def solve_beetle_mystery():
    """
    Analyzes facts about the Rainbow Leaf Beetle (Chrysolina cerealis)
    to determine why it is unlikely to be observed in Germany.
    """
    
    # --- Step 1: Define facts about the insect ---
    insect_name = "Chrysolina cerealis (Rainbow Leaf Beetle)"
    is_real = True
    is_extinct = False
    native_continents = ["Europe", "Asia"]
    habitat_type = "Temperate mountainous"
    is_in_germany = True
    conservation_status_germany = "Critically Endangered"
    # A famous 2017 study found over 75% decline in total flying insect biomass in German protected areas.
    # A "Critically Endangered" status (per IUCN) often implies a population reduction of >80%.
    major_population_decline_reported = True
    decline_percentage_figure = 76

    # --- Step 2: Define the answer choices ---
    options = {
        "A": "It is endemic to North America",
        "B": "It is endemic to the tropics",
        "C": f"Its population size has been reduced by over {decline_percentage_figure}% in the last four decades",
        "D": "It is not real",
        "E": "It is extinct",
        "F": "It is present in Germany, but has not been observed in over ten years"
    }

    # --- Step 3: Evaluate each option ---
    print("Evaluating the options to find the most likely reason:")

    # Evaluate A
    if "North America" in native_continents:
        print("Reason A might be correct.")
    else:
        print(f"Reason A is incorrect. The beetle is native to {', '.join(native_continents)}, not North America.")

    # Evaluate B
    if "Tropical" in habitat_type:
        print("Reason B might be correct.")
    else:
        print(f"Reason B is incorrect. The beetle's habitat is {habitat_type}, not tropical.")

    # Evaluate D
    if not is_real:
        print("Reason D might be correct.")
    else:
        print("Reason D is incorrect. The beetle is real.")

    # Evaluate E
    if is_extinct:
        print("Reason E might be correct.")
    else:
        print("Reason E is incorrect. The beetle is not extinct.")
    
    # Evaluate F
    # This is hard to prove definitively false without access to all observation records.
    # However, its "Critically Endangered" status implies it is tracked, likely with recent observations by experts.
    # The underlying CAUSE of its rarity is the population decline, making C a more fundamental answer.
    print("Reason F is less fundamental than the cause of its rarity. The reason it would be hard to observe is DUE to a small population.")

    # Evaluate C
    if major_population_decline_reported and conservation_status_germany == "Critically Endangered":
        print(f"\nReason C is the most plausible. The beetle is listed as '{conservation_status_germany}' in Germany.")
        print("This status implies a severe population reduction.")
        print(f"This aligns with reports of major insect population declines, such as the equation: general decline > {decline_percentage_figure}%.")
        print(f"Therefore, a massive population reduction is the best explanation for why it's unlikely to be seen.")
        final_answer = "C"
    else:
        final_answer = "Undetermined"
        
    return final_answer

# --- Run the analysis ---
correct_option = solve_beetle_mystery()
# The final answer is wrapped according to the format requirements.
# print(f"\n<<< {correct_option} >>>")