def final_verification(combo):
    # Helper function to check if a letter is between two other letters
    def is_between(letter, early, late):
        return ord(early) < ord(letter) < ord(late)

    # We know from condition 7,0,'K','G' that K is too early in the alphabet
    # This is a crucial condition we might have missed earlier
    if combo[3] == 'K':
        # From condition 14: 70KG - both letters are incorrect and too early
        # This explicitly states K is too early in the alphabet
        return False
    
    # Additional verification for M
    if combo[3] == 'M':
        # Check if M satisfies all position requirements:
        
        # Must be after early letters (E,F,H,I,B,K,G)
        early_letters = ['E','F','H','I','B','K','G']
        for letter in early_letters:
            if not ord('M') > ord(letter):
                return False
                
        # Must be before late letters (P,Z,V,W,X,N,Y)
        late_letters = ['P','Z','V','W','X','N','Y']
        for letter in late_letters:
            if not ord('M') < ord(letter):
                return False
        
        # M is never explicitly mentioned as incorrect in any guess
        # M is in the middle range of the alphabet
        # M satisfies all relative position requirements
        return True
    
    return False

# Test both remaining candidates
candidates = [[4, 5, 'L', 'K'], [4, 5, 'L', 'M']]
print("Final verification:")
for candidate in candidates:
    print(f"{candidate}: {'Valid' if final_verification(candidate) else 'Invalid'}")