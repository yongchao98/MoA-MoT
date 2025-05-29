def check_peng(hand, new_card):
    return hand.count(new_card) >= 2

def check_chi(hand, new_card):
    sorted_hand = sorted(hand + new_card)
    for i in range(len(sorted_hand) - 2):
        if ord(sorted_hand[i+1]) == ord(sorted_hand[i]) + 1 and ord(sorted_hand[i+2]) == ord(sorted_hand[i]) + 2:
            return True
    return False

def determine_result(hand, new_card):
    if check_peng(hand, new_card):
        return 1  # Peng
    elif check_chi(hand, new_card):
        return 2  # Chi
    else:
        return 0  # Pass

# Initial hand
hand = list("AVRNHDROSGQYP")

# Round 1
hand.remove('H')
hand.append('I')
result1 = determine_result(hand, 'I')

# Round 2
hand.remove('V')
hand.append('T')
result2 = determine_result(hand, 'T')

# Round 3
hand.remove('Y')
hand.append('R')
result3 = determine_result(hand, 'R')

# Round 4
hand.remove('P')
hand.append('S')
result4 = determine_result(hand, 'S')

# Round 5
hand.append('M')
result5 = determine_result(hand, 'M')

print(result5)