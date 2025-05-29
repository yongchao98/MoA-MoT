def calculate_score(s):
    score = 0
    for i in range(len(s) - 1):
        if s[i] == s[i + 1]:
            if s[i] == 'A':
                score += 1
            elif s[i] == 'B':
                score += 2
            elif s[i] == 'C':
                score += 3
            elif s[i] == 'D':
                score += 4
            elif s[i] == 'E':
                score += 5
            elif s[i] == 'F':
                score += 6
    return score

string = "FDFDADAEBCAEAFCCAEDFDEDAA"
total_score = calculate_score(string)
print(total_score)