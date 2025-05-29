from collections import defaultdict

def most_active(bio_data):
    active_counts = defaultdict(int)
    
    for name, start, end in bio_data:
        for year in range(start, end + 1):
            active_counts[year] += 1
    
    max_count = max(active_counts.values())
    max_years = [year for year, count in active_counts.items() if count == max_count]
    
    start_year = min(max_years)
    end_year = max(max_years)
    
    return (start_year, end_year)

# Test the function with a constructed input
bio_data = [
    ["Author A", 1990, 1991],
    ["Author B", 1991, 1992],
    ["Author C", 1991, 1991]
]

result = most_active(bio_data)
print(result)