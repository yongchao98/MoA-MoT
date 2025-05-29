from collections import Counter

def can_form_palindrome(letters):
    count = Counter(letters)
    odd_count = sum(1 for c in count if count[c] % 2 != 0)
    return odd_count <= 1

letters = ['o', 's', 'o', 'm', 'm', 's']
if can_form_palindrome(letters):
    palindrome = "osmmso"
    print(palindrome)
else:
    print("Cannot form a palindrome")