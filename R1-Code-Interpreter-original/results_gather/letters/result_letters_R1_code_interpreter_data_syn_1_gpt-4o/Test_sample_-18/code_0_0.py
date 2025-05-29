word = 'nmfghmzgvjinzgmijkdxs'
g_positions = [i + 1 for i, char in enumerate(word) if char == 'g']
g_count = len(g_positions)
print(f"Count: {g_count}, Positions: {g_positions}")