from collections import Counter

def form_palindrome(letters):
    count = Counter(letters)
    half_palindrome = []
    middle = []

    for letter, cnt in count.items():
        if cnt % 2 != 0:
            middle.append(letter)
        half_palindrome.append(letter * (cnt // 2))

    # If more than one character has an odd count, a palindrome is not possible
    if len(middle) > 1:
        return None

    # Form the palindrome
    half = ''.join(half_palindrome)
    return half + ''.join(middle) + half[::-1]

letters = ['w', 'w', 'f', 'f']
palindrome = form_palindrome(letters)
print(palindrome)