# Import the regular expression module
import re

# Constructing an input where no address meets the criteria for part1 or part2
lines = """
abcd[efgh]ijkl
mnop[qrst]uvwx
yzab[cdef]ghij
klmn[opqr]stuv
"""

# None of these lines have an ABBA in the outer block.
# None of these lines have an ABA in the outer block with a corresponding BAB in the inner block.

def splitBlocks(line):
    inners = []
    outers = []
    for idx, block in enumerate(re.split(r'[\[\]]', line.strip())):
        if idx % 2 == 0:
            outers.append(block)
        else:
            inners.append(block)
    return [inners, outers]

def containsABBA(txt):
    for idx in range(0, len(txt)-3):
        a, b, c, d = txt[idx:idx+4]
        if a == d and b == c and a != b:
            return True
    return False

def listABAs(txt):
    for idx in range(0, len(txt)-2):
        a, b, c = txt[idx:idx+3]
        if a == c and a != b:
            yield txt[idx:idx+3]

def invertABA(aba):
    return aba[1] + aba[0] + aba[1]

def main_solution(lines):
    lines = lines.split('\n')
    count_part1 = 0
    for line in lines:
        inners, outers = splitBlocks(line.strip())
        if any(filter(containsABBA, outers)) and not any(filter(containsABBA, inners)):
            count_part1 += 1
    
    count_part2 = 0
    for line in lines:
        inners, outers = splitBlocks(line.strip())
        for bab in [invertABA(aba) for outer in outers for aba in listABAs(outer)]:
            if any(filter(lambda inner: bab in inner, inners)):
                count_part2 += 1
                break
    
    return {"part1": count_part1, "part2": count_part2}

# Execute the solution with the constructed input
result = main_solution(lines)
print(result)