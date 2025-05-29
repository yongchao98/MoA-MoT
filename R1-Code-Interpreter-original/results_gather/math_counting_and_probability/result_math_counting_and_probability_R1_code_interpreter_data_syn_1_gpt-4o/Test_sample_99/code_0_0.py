# Calculate the expected value of a fair octahedral die
faces = list(range(1, 9))
probability = 1 / 8
expected_value = sum(face * probability for face in faces)
print(expected_value)