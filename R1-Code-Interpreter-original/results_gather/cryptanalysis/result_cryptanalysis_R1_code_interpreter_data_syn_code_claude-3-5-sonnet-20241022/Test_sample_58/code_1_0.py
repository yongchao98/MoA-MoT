def check_letter_sequence(first_letter):
    # From #3: G and E are incorrect and too early - means our letter should be after G
    # From #5: L and E are incorrect and too early - means our letter should be after L
    # From #7: F is incorrect and too early - means our letter should be after F
    # Therefore, our first letter should be after L in the alphabet
    
    # From the remaining letters [M,O,T,U,V,W,X,Y], we need one that:
    # 1. Is after L (already satisfied)
    # 2. Makes sense with all the feedback about "too early in alphabet"
    # 3. Works with Z being confirmed correct in second position
    
    # Convert letter to ASCII value for comparison
    letter_val = ord(first_letter)
    
    # Must be after L (ord('L') = 76)
    if letter_val <= ord('L'):
        return False
        
    # Test against all the guesses where letters were "too early"
    if letter_val <= ord('G') or letter_val <= ord('L'):
        return False
    
    return True

candidates = [['5', '7', l1, 'Z'] for l1 in 'ABDHIJKMOTUVWXY']
valid_ones = [cand for cand in candidates if check_letter_sequence(cand[2])]
print(valid_ones)