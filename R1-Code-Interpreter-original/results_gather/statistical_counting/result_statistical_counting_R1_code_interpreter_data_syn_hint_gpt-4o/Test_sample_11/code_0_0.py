def calculate_score(s):
    total_score = 0
    i = 0
    while i < len(s) - 1:
        if s[i] == s[i + 1]:
            if s[i] == 'A':
                total_score += 1
            elif s[i] == 'B':
                total_score += 2
            elif s[i] == 'C':
                total_score += 3
            elif s[i] == 'D':
                total_score += 4
            elif s[i] == 'E':
                total_score += 5
            # Skip the next character as it is part of the consecutive sequence
            i += 1
        i += 1
    return total_score

string = "ABCEEBEBCEECEDDBEBCA"
score = calculate_score(string)
print(score)