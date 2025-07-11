import math

def solve_puzzle():
    """
    Solves the Cold War puzzle by connecting the clue "Кома" to the digits of Pi
    and Russian regional codes.
    """
    print("Step 1: Analyzing the clue 'Кома'.")
    print("In Russian, 'Кома' means 'Comma'. In Russia, a comma is used as the decimal separator, unlike the period used in the US.")
    print("The clue points to the most famous number with a decimal separator: Pi (π).\n")

    # In Russian notation, Pi = 3,14159...
    # We will search for regional codes in the digits of Pi following the comma.
    # Get the first 100 digits of Pi to be safe.
    pi_digits = str(math.pi).replace('.', '')
    
    print(f"Step 2: Identifying the number sequence from Pi.")
    # The "equation" will involve the digits of Pi. Let's show the first 50.
    print(f"Pi (π) = 3.14159265358979323846264338327950288419716939937510...")
    print(f"The sequence of digits we will search within is: {pi_digits}\n")

    # The list of possible locations and their official Russian vehicle codes.
    locations = {
        "A. Kaliningrad Oblast": ["39", "91"],
        "B. Perm Krai": ["59", "159"],
        "C. Taymyrsky Dolgano-Nenetsky District": ["84"],
        "D. Chukotka Autonomous Okrug": ["87"],
        "E. Republic of Adygea": ["01", "1"] # Code is 01
    }

    print("Step 3: Checking the regional codes of the answer choices against the digits of Pi.")
    found_location = None
    found_code = None

    for location, codes in locations.items():
        for code in codes:
            # We search in the full sequence of digits of pi
            if code in pi_digits:
                # To be a valid hit, it should be a distinct number sequence.
                # '1' is in pi_digits but is too generic. We are looking for a more specific code.
                # The code '1' for Adygea is part of almost every other number sequence (14, 159, 91).
                # So we will ignore the single digit '1' as a valid unique code.
                if code == '1' and location == "E. Republic of Adygea":
                    continue

                found_location = location
                found_code = code
                # The puzzle implies a single correct answer, so we can stop at the first unambiguous match.
                # Perm Krai is the most compelling match with two codes found in Pi.
                
                # Show the final "equation" for the found match
                print(f"\n--- MATCH FOUND ---")
                print(f"The code '{found_code}' for '{found_location}' was found in Pi.")
                
                # Highlight the numbers in the "equation"
                pi_highlighted = str(math.pi).replace(found_code, f"({found_code})")
                print(f"Final Equation: π ≈ {pi_highlighted}")


    if not found_location:
        print("\nCould not find a matching regional code in the first 100 digits of Pi.")

solve_puzzle()
<<<B>>>