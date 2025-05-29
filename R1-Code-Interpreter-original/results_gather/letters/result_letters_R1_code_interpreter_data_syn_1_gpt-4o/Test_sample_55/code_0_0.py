word = "pseudoparthenogenesis"
count_g = word.count('g')
positions = [i + 1 for i, letter in enumerate(word) if letter == 'g']
print(f"Count: {count_g}, Positions: {positions}")