s = "bbbcccbcb"
print(f"Starting string: {s}")
print(f"Number of b's: {s.count('b')}")
if s.count('b') > 0 and s.count('b') % 2 == 0:
    final = s + "ab"
    print(f"After adding 'ab': {final}")
else:
    final = s
    print("No change needed")
print(f"\nFinal string: {final}")