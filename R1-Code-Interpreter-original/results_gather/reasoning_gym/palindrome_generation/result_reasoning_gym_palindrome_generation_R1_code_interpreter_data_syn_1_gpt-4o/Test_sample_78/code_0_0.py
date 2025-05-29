from collections import Counter

# List of letters
letters = ['u', 'p', 'u', 's', 'k', 'k', 's', 'u', 'p']

# Count the frequency of each letter
frequency = Counter(letters)

# Check the number of odd frequency letters
odd_count = sum(1 for count in frequency.values() if count % 2 != 0)

# Determine if a palindrome can be formed
can_form_palindrome = odd_count <= 1

# If a palindrome can be formed, construct it
if can_form_palindrome:
    half_palindrome = []
    middle_char = ''
    
    for char, count in frequency.items():
        if count % 2 == 0:
            half_palindrome.append(char * (count // 2))
        else:
            half_palindrome.append(char * (count // 2))
            middle_char = char
    
    # Form the palindrome
    first_half = ''.join(half_palindrome)
    palindrome = first_half + middle_char + first_half[::-1]
    print(palindrome)
else:
    print("Cannot form a palindrome")