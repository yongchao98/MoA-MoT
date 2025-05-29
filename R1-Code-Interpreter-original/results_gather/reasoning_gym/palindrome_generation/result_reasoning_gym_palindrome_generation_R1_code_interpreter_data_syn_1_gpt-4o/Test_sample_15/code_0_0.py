from collections import Counter

def form_palindrome(letters):
    count = Counter(letters)
    half_palindrome = []
    middle = []

    for letter, cnt in count.items():
        if cnt % 2 == 1:
            middle.append(letter)
        half_palindrome.extend([letter] * (cnt // 2))

    # Sort to ensure consistent output
    half_palindrome.sort()
    middle.sort()

    # Form the palindrome
    return ''.join(half_palindrome + middle + half_palindrome[::-1])

letters = ['o', 'g', 'g', 'b', 'b', 'o', 'o', 'o']
palindrome = form_palindrome(letters)
print(palindrome)