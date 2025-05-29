def check_numbers(n1, n2):
    # Check number conditions
    if n1 < 0 or n2 < 0 or n1 > 9 or n2 > 9:
        return False
    if n1 == n2:  # no repeating
        return False
    if n1 <= 0 and n2 <= 3:  # condition 1
        return False
    if n1 >= 7 and n2 >= 6:  # condition 2
        return False
    if (n1 == 8 and n2 == 3) or (n1 == 3 and n2 == 8):  # condition 3
        return False
    if n1 >= 8 or n2 >= 8:  # condition 5
        return False
    return True

def check_letters(l1, l2):
    # Check letter conditions
    if l1 == l2:  # no repeating
        return False
    if l1 in 'OQ' and l2 in 'OQ':  # condition 1
        return False
    if l1 in 'EC' and l2 in 'EC':  # condition 2
        return False
    return True

def find_solution():
    numbers = range(10)
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    # Try all possible combinations
    for n1 in numbers:
        for n2 in numbers:
            if not check_numbers(n1, n2):
                continue
            
            # If we get here, the numbers satisfy all conditions
            for l1 in letters:
                for l2 in letters:
                    if not check_letters(l1, l2):
                        continue
                    
                    # Check condition 4
                    if ((n1 == 4 or n2 == 4) and  # one number is 4
                        (l1 in 'OF' and l2 in 'OF')):  # one letter is O or F
                        
                        # Check condition 3
                        if ((l1 == 'F' and l2 != 'O') or 
                            (l1 != 'O' and l2 == 'F')):
                            
                            # Check condition 5
                            if ((l1 == 'L' and l2 != 'K') or 
                                (l1 != 'K' and l2 == 'L')):
                                
                                print([str(n1), str(n2), l1, l2])

find_solution()