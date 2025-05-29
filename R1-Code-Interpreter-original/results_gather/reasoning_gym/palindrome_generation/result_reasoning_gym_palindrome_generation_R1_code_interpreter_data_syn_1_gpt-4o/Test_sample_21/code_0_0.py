from collections import Counter

def can_form_palindrome(letters):
    count = Counter(letters)
    odd_count = sum(1 for c in count.values() if c % 2 != 0)
    return odd_count <= 1

letters = ['s', 'y', 'y', 'l', 'l', 's']
if can_form_palindrome(letters):
    palindrome = "sylls"
    print(palindrome)
else:
    print("Cannot form a palindrome")