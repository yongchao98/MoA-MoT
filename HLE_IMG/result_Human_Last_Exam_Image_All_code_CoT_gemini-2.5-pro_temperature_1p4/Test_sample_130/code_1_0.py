def identify_element_from_spectrum(observed_lines, elements_db, tolerance=5.0):
    """
    Identifies an element by matching observed spectral lines against a database.

    Args:
        observed_lines (list): A list of wavelengths (in nm) for the observed lines.
        elements_db (dict): A dictionary with element names as keys and lists of their
                              prominent spectral lines (in nm) as values.
        tolerance (float): The margin of error (in nm) for a match.

    Returns:
        str: The name of the element with the most matches.
    """
    match_counts = {element: 0 for element in elements_db}
    
    print("Starting spectral analysis...")
    print(f"Observed Lines (nm): {observed_lines}\n")

    for element, known_lines in elements_db.items():
        print(f"--- Checking against {element} ---")
        matches = []
        # To avoid double-counting, we track which known lines have been matched.
        matched_known_lines = [False] * len(known_lines)
        
        for obs_line in observed_lines:
            for i, known_line in enumerate(known_lines):
                if abs(obs_line - known_line) <= tolerance and not matched_known_lines[i]:
                    match_counts[element] += 1
                    matches.append(f"Observed {obs_line}nm -> {element} {known_line}nm")
                    matched_known_lines[i] = True
                    # A single observed line should match at most one known line.
                    break
        
        if matches:
            for match in matches:
                print(f"  [MATCH] {match}")
        print(f"Total matches for {element}: {match_counts[element]}\n")

    # Find the element with the highest match count
    if not any(match_counts.values()):
        return "No element identified"
        
    best_match_element = max(match_counts, key=match_counts.get)
    return best_match_element

# --- Data ---
# A simplified database of prominent visible spectral lines (in nm).
# Iron (Fe) has a very rich spectrum; this is a small subset.
SPECTRAL_LINES_DB = {
    "Hydrogen": [410.2, 434.1, 486.1, 656.3],
    "Helium": [447.1, 471.3, 492.2, 501.6, 587.6, 667.8],
    "Iron": [438.4, 440.5, 495.7, 527.0, 537.1, 542.9, 619.2, 639.4, 640.0],
    "Neon": [585.2, 588.2, 603.0, 607.4, 614.3, 621.7, 633.4, 640.2, 650.6]
}

# Wavelengths estimated from the prominent lines in the provided image.
# The image shows many lines, especially clustered in blue and green.
OBSERVED_PROMINENT_LINES = [
    440,  # Prominent blue line
    496,  # Cyan/Green line
    527,  # Very strong green line cluster
    620,  # Orange-red line
    640   # Strong red line cluster
]


# --- Execution ---
identified_element = identify_element_from_spectrum(OBSERVED_PROMINENT_LINES, SPECTRAL_LINES_DB)

print("="*30)
print(f"Conclusion: The element with the most matches is {identified_element}.")
print("="*30)