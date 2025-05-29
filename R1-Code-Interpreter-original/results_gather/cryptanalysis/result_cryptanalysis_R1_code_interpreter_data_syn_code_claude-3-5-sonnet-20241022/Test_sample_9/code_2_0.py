def analyze_numbers():
    # All numbers that appeared in guesses
    tested_numbers = set(['6', '4', '5', '1', '2', '9', '7'])
    
    # Possible first digits
    candidates = ['0', '3', '8']
    
    print("Analysis of candidate first digits:")
    for num in candidates:
        # If a number was never tested in any guess, we can't be certain about it
        if num not in tested_numbers:
            print(f"{num}: Cannot be verified as we have no information about it from any guess")
        else:
            print(f"{num}: Was tested in guesses")
            
    # Check which numbers were proven wrong
    wrong_numbers = set(['6', '4', '5', '1', '2', '9'])
    print("\nNumbers proven wrong:", wrong_numbers)
    
    # Numbers that were tested and not proven wrong
    valid_numbers = tested_numbers - wrong_numbers
    print("Numbers tested and not proven wrong:", valid_numbers)
    
    # Numbers never tested
    untested_numbers = set(['0', '3', '8']) - tested_numbers
    print("Numbers never tested:", untested_numbers)

analyze_numbers()