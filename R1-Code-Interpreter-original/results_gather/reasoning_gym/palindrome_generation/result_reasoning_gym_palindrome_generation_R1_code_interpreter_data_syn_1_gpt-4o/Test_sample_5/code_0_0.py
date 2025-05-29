from collections import Counter

def can_form_palindrome(letters):
    count = Counter(letters)
    odd_count = sum(1 for c in count.values() if c % 2 != 0)
    return odd_count <= 1

letters = ['y', 'r', 'r', 'y']
if can_form_palindrome(letters):
    palindrome = "yrry"
    print(palindrome)
else:
    print("Cannot form a palindrome")